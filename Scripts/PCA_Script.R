# A script to create PCAs for 24 hour biogeochem data for Cabral and Varari
# Edited on 6/10/2023
# Created by Nyssa Silbiger 

####################################

# load libraries #################
library(tidyverse)
library(ggfortify)
library(lubridate)
library(ggforce)
library(viridis)
library(patchwork)
library(here)
library(wesanderson)
library(broom)
library(ggtext)
library(glue)
library(htmlTable)

# load the 24 hour chemistry data #####################
Data<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/August2021/Allbiogeochemdata_QC2.csv") %>% mutate(Season = "Dry")

Data_march<-read_csv("https://raw.githubusercontent.com/njsilbiger/MooreaSGD_site-selection/main/Data/March2022/Allbiogeochemdata_QC_march_fdom2.csv")%>% mutate(Season = "Wet") %>%
  mutate(DateTime = mdy_hms(paste(Date,as.character(Time))))


# Load the Seep Data
VarariSeep<-read_csv(here("Data","Varari","AllVarariSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site, Season = Season_seep)

CabralSeep<-read_csv(here("Data","Cabral","AllCabralSeepData.csv")) %>%
  rename_with(.cols = TempInSitu:PAR_calc,function(x){paste0(x,"_seep")}) %>% # rename columns to say seep at the end
  rename(DateTime = date, Location = Site, Season = Season_seep)

# bind them together
SeepAll<-bind_rows(VarariSeep, CabralSeep)

# Bind with the chem data
Data<-Data %>%
  bind_rows(Data_march)%>%
  left_join(SeepAll) %>%
   mutate(NewDay_Night = case_when( 
    TimeBlock == "Evening" ~ "Day", # change day and night delineations... dusk and dawn are now night and day 
    TimeBlock == "Morning" ~ "Night" ,
    TimeBlock == "Night" ~ "Night",
    TimeBlock == "Afternoon" ~ "Day"
    # Tide == "Low" & Day_Night == "Day" ~ "Night", # change day and night delineations... dusk and dawn are now night and day 
    # Tide == "Low" & Day_Night == "Night" ~ "Day",
    # Tide == "High" & Day_Night == "Day" ~ "Day",
    # Tide == "High" & Day_Night == "Night" ~ "Night"
  ),
  Day_Night = NewDay_Night) %>% # replace it
  select(!NewDay_Night)



### Varari #####
## There seems to be a contaminated nutrient sample for V2 Low tide on the 8/8/2021.  Remove this point
remove<-Data %>% filter(CowTagID=="V2", Tide == "Low",Day_Night=="Night", Date =="8/8/2021") 
removea<-Data %>% filter(CowTagID == "V17", Tide == "High", Day_Night == "Night", Season == "Wet")# Ammonia is an outlier

## also filter out all the data from the first low tide where the water level was super high
remove2<-Data %>% filter(Tide =="Low", Day_Night=="Night", Date == "8/6/2021")
remove_varari<-bind_rows(remove, removea, remove2) # bring the "bad data" into one tibble
 
# fdom or silicate data off for these
remove3<-Data %>% filter(CowTagID== "C4", Tide =="Low", Day_Night=="Day", Date == "8/9/2021")
remove4<-Data %>% filter(CowTagID== "C2", Tide =="High", Day_Night=="Night", Date == "3/31/2022")
remove5<-Data %>% filter(CowTagID== "C17", Tide =="Low", Day_Night=="Day", Date == "3/30/2022")

remove_cabral<- bind_rows(remove3, remove4,remove5) # bring the "bad data" into one tibble

## create log transformed data for everything except for salinity, temperature, TA, and pH
Datalog<-Data %>%
  mutate(Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tyrosine_Like+Tryptophan_Like)%>%
  mutate_at(vars(Phosphate_umolL:Lignin_Like, Humics, Proteinaceous), .funs = function(x){log(x+0.001)})

remove_vararilog<-remove_varari %>% mutate_at(vars(Phosphate_umolL:Lignin_Like), .funs = function(x){log(x+0.001)})
remove_cabrallog<-remove_cabral %>% mutate_at(vars(Phosphate_umolL:Lignin_Like), .funs = function(x){log(x+0.001)})


# extract the params for the PCA
V_pca_Data_both<-Datalog %>%
   anti_join(remove_vararilog)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  select(Season, Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, Humics,Proteinaceous, TA)%>%
  group_by(Season)%>%
  mutate_all(.funs = scale)%>% # scale by season
  drop_na() %>%
  ungroup()%>%
  select(-Season)

# Run the PCA
pca_V_both <- prcomp(V_pca_Data_both)


# calculate percent explained by each PC
perc.explained_both<-round(100*pca_V_both$sdev/sum(pca_V_both$sdev),1)

# Extract the scores and loadings
PC_scores_both <-as_tibble(pca_V_both$x[,1:2])


PC_loadings_both<-as_tibble(pca_V_both$rotation) %>%
  bind_cols(labels = rownames(pca_V_both$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like", "MarineHumic_Like", "Humics","Proteinaceous")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        #   labels == "Lignin_Like" ~"Lignin Like",
                        #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity",
                        labels == "Humics" ~ "Humic fDOM",
                        labels == "Proteinaceous" ~"Proteinaceous fDOM"))

# Put it with all the original data

V_pca_Data_all_both<-Data %>%
  anti_join(remove_varari)%>%
  #anti_join(remove2)%>%
  filter(Location == "Varari", Plate_Seep=="Plate") %>%
  drop_na(Salinity,pH,Phosphate_umolL:Lignin_Like )%>%
  bind_cols(PC_scores_both)

# scores plot
p1_both<-V_pca_Data_all_both %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = TimeBlock))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#D64550","#EA9E8D"))+
  scale_fill_manual(values = c("#D64550","#EA9E8D"))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  ggforce::geom_mark_ellipse(
    aes(#fill = Tide, 
        label = paste(TimeBlock, Tide), color =Tide), 
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap=0, tol = 0.05)+
  geom_point(size = 2) +
  labs(
      # x = paste0("PC1 ","(",perc.explained_both[1],"%)"),
       x = "",
       y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        strip.text = element_blank())+
  facet_wrap(~Season, nrow=2)

# loadings plot 
p2_both<-PC_loadings_both %>%
  ggplot(aes(x=PC1, y=PC2, label=nicenames, color = groupings))+
  geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA) +
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_segment(data = PC_loadings_both, aes(x=0,y=0,xend=PC1*10,yend=PC2*10),size = 1.2,
               arrow=arrow(length=unit(0.1,"cm")))+
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(color ="",
       y = "",
       x = paste0("PC1 ","(",perc.explained_both[1],"%)"))+
       #y = paste0("PC2 ","(",perc.explained_both[2],"%)"))+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.75, 0.75),
        legend.position = "none",
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16))


## bring it together
VarariPCA_botha<-p1_both/plot_spacer()/p2_both+plot_layout(heights = c(2,-0.1, 1))

# all PCA
VarariPCA_both<-(p1_both)/(p2_both)+ 
  patchwork::plot_annotation(#"Cabral Plates", 
    theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", 
                                            hjust = 0.5, 
                                            margin = margin(t = 10, b = 20, 
                                                            unit = "pt"))))

ggsave(plot = VarariPCA_both, filename = here("Output","VarariPCA_both.pdf"), width = 15, height = 15)
#### Cabral #####

#Extract the cabral data

C_pca_Data_both<-Datalog %>%
  anti_join(remove_cabrallog)%>%
  filter(Location == "Cabral", Plate_Seep=="Plate") %>%
  select(Season, Salinity,pH,Phosphate_umolL:NN_umolL,Ammonia_umolL, Humics,Proteinaceous, TA)%>%
  group_by(Season)%>%
  mutate_all(.funs = scale)%>% # scale by season
  drop_na() %>%
  ungroup()%>%
  select(-Season)

pca_C_both <- prcomp(C_pca_Data_both)

# calculate percent explained by each PC
perc.explainedC_both<-round(100*pca_C_both$sdev/sum(pca_C_both$sdev),1)


# Extract the scores and loadings
PC_scoresC_both <-as_tibble(pca_C_both$x[,1:2])

PC_loadingsC_both<-as_tibble(pca_C_both$rotation) %>%
  bind_cols(labels = rownames(pca_C_both$rotation))%>%
  mutate(groupings = case_when( # add groupings
    labels %in% c("Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL")~ "Nutrient Chemistry",
    labels == "Salinity" ~ "Salinity",
    labels %in% c("pH","TA") ~ "Carbonate Chemistry",
    labels %in% c("HIX","Tryptophan_Like","Tyrosine_Like","VisibleHumidic_Like", "MarineHumic_Like", "Humics","Proteinaceous")~"fDOM"
  ),
  nicenames = case_when(labels == "TempInSitu_seep" ~ "Temperature",
                        labels == "pH" ~ "pH<sub>T</sub>",
                        #   labels == "Lignin_Like" ~"Lignin Like",
                        #    labels == "M_C" ~ "M:C",
                        labels == "Tyrosine_Like" ~ "Tyrosine Like",
                        labels == "Tryptophan_Like" ~ "Tryptophan Like",
                        labels == "HIX"~"HIX",
                        labels == "MarineHumic_Like" ~ "Marine Humic Like",
                        labels == "VisibleHumidic_Like" ~ "Visible Humic Like",
                        labels == "Ammonia_umolL" ~ "Ammonium",
                        labels == "TA" ~ "Total Alkalinity",
                        labels == "Phosphate_umolL" ~ "Phosphate",
                        labels == "NN_umolL" ~ "Nitrate+Nitrite",
                        labels == "Silicate_umolL" ~ "Silicate",
                        labels == "Salinity" ~"Salinity",
                        labels == "Humics" ~ "Humic fDOM",
                        labels == "Proteinaceous" ~"Proteinaceous fDOM"))

# Put it with all the original data

C_pca_Data_all_both<-Data %>%
  anti_join(remove_cabral)%>%
  select(!Jamie_Plate_ID)%>% # Jamie's plates are all NA here
  filter(Location == "Cabral", Plate_Seep=="Plate") %>%
  drop_na(Salinity,pH,Phosphate_umolL:NN_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like,TA) %>%
  bind_cols(PC_scoresC_both)  

# scores plot

p1c_both<-C_pca_Data_all_both %>%
  ggplot(aes(x = PC1, y = PC2, color = Tide, shape = TimeBlock))+
  geom_hline(aes(yintercept = 0), lty = 2)+
  geom_vline(aes(xintercept = 0), lty = 2)+
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  scale_shape_manual(values = c(1, 22,15,16))+
  scale_colour_manual(values = c("#16697A","#82C0CC"))+
  scale_fill_manual(values = c("#16697A","#82C0CC"))+
  labs(
     #  x = paste0("PC1 ","(",perc.explainedC_both[1],"%)"),
       y = paste0("PC2 ","(",perc.explainedC_both[2],"%)"),
       x = "")+
  ggforce::geom_mark_ellipse(
    aes(#fill = Tide, 
        label = paste(TimeBlock, Tide), color = Tide), 
    alpha = .35, show.legend = FALSE,  label.buffer = unit(1, "mm"), con.cap = 0)+
  geom_point(size = 2) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.background = element_blank(),
        strip.text = element_blank())+
  facet_wrap(~Season, nrow=2)



# loadings plot 

p2c_both<-PC_loadingsC_both %>%
  ggplot(aes(x=PC1, y=PC2, label=nicenames, color = groupings))+
  geom_hline(aes(yintercept = 0), lty = 2)+
  geom_vline(aes(xintercept = 0), lty = 2)+
  geom_segment(aes(x=0,y=0,xend=PC1*10,yend=PC2*10),
               arrow=arrow(length=unit(0.1,"cm")), size = 1.2)+
   geom_richtext(aes(x = PC1*10+0.1, y = PC2*10+.1 ), show.legend = FALSE, size = 5, fill=NA, label.colour = NA)+
  scale_color_manual(values = wes_palette("Darjeeling1"))+
   coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(color = "",
       x = paste0("PC1 ","(",perc.explainedC_both[1],"%)"),
       #y = paste0("PC2 ","(",perc.explainedC_both[2],"%)"),
       y = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.75, 0.80),
        legend.text = element_markdown(size = 16),
        legend.key.size = unit(1, 'cm'),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.background = element_blank())


CabralPCA_botha<-p1c_both/plot_spacer()/p2c_both+plot_layout(heights = c(2,-0.1, 1))


VarariPCA_botha|CabralPCA_botha
ggsave(here("Output","BothSites_PCA.pdf"), height = 15, width = 15)


# all PCA
CabralPCA_both<-(p1c_both)/(p2c_both)+ 
  patchwork::plot_annotation(#"Cabral Plates", 
    theme = theme(plot.title = element_text(size = rel(1.5), face = "bold", 
                                            hjust = 0.5, 
                                            margin = margin(t = 10, b = 20, 
                                                            unit = "pt"))))

ggsave(plot = CabralPCA_both, filename = here("Output","CabralPCA_both.pdf"), width = 15, height = 15)


## Calculate the ranges for each parameter
ranges<-Data %>%
  filter(Plate_Seep=="Seep"| Plate_Seep == "Spring") %>% # just do the seep data
  filter(CowTagID != "Varari_Well") %>% # remove the well sample
  mutate(Humics = MarineHumic_Like+VisibleHumidic_Like,
         Proteinaceous = Tryptophan_Like+Tyrosine_Like)%>%
  select(Location,Season, Salinity, TA,pH,Phosphate_umolL:Ammonia_umolL, VisibleHumidic_Like, Tyrosine_Like, Tryptophan_Like, HIX, TempInSitu_seep, MarineHumic_Like, Lignin_Like, Humics, Proteinaceous, M_C)%>%
 # select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep) %>%
  pivot_longer(cols = c(Salinity, TA: M_C)) %>%
  drop_na()%>%
  group_by(Location,name, Season)%>%
  summarise(min = round(min(value),2),
            max = round(max(value),2)) %>%
    mutate(unit = case_when(name == "TempInSitu_seep" ~ " &deg;C", # add units
                            name %in% c("Ammonia_umolL","NN_umolL","Silicate_umolL","Phosphate_umolL")~" &mu;mol L<sup>-1</sup>",
                            name %in% c("M_C","HIX")~ " Ratio",
                            name %in% c("Salinity")~ " psu",
                            name %in% c("pH")~ " pH total scale",
                            name == "TA"~ " &mu;mol kg<sup>-1</sup>",
                            name %in% c("Lignin_Like","Tyrosine_Like","Tryptophan_Like","MarineHumic_Like","VisibleHumidic_Like","Humics","Proteinaceous")~" Raman Units"),
      range = paste0("[",min," - ",max,unit,"]"),
      nicenames = case_when(name == "TempInSitu_seep" ~ "Temperature",
                            name == "pH" ~ "pH<sub>T</sub>",
                            name == "Lignin_Like" ~"Lignin Like",
                         #   name == "M_C" ~ "M:C",
                            name == "Tyrosine_Like" ~ "Tyrosine-Like",
                            name == "Tryptophan_Like" ~ "Tryptophan-Like",
                            name == "HIX"~"HIX",
                            name == "MarineHumic_Like" ~ "Marine Humic-Like",
                            name == "VisibleHumidic_Like" ~ "Visible Humic-Like",
                            name == "Ammonia_umolL" ~ "Ammonium",
                            name == "TA" ~ "Total Alkalinity",
                            name == "Phosphate_umolL" ~ "Phosphate",
                            name == "NN_umolL" ~ "Nitrate+Nitrite",
                            name == "Silicate_umolL" ~ "Silicate",
                            name == "Salinity" ~"Salinity",
                            name == "Humics"~ "Humics fDOM",
                            name == "Proteinaceous"~"Proteinaceous fDOM"
                            
        
      ))
  


# Make a correlation plot of everything versus Salinity

cor_fun <- function(data) cor.test(data$value, data$Salinity, method = "pearson")%>% tidy()

cortest<-Datalog %>%
  filter(Plate_Seep=="Seep" | Plate_Seep == "Spring" ) %>% # just do the seep data
  filter(CowTagID != "Varari_Well") %>% # remove the well sample
  select(Location, Season, Salinity, TA: Lignin_Like, Humics, Proteinaceous, TempInSitu_seep, Depth_seep)%>%
  pivot_longer(cols = c(TA: Lignin_Like, Humics, Proteinaceous, TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(Location, name, Season)%>%
  nest() %>%
  mutate(model = map(data, cor_fun)) %>%
  select(Location, name,Season, model) %>%
  unnest(model)%>% # calculate correlation coefficient
  mutate(sig = ifelse(p.value<0.05, 1,0 ))# add a 1 if significant correlation (to 0.1 pvalue)


# Plot it
cortest %>%
  left_join(ranges) %>% # join with the ranges
  mutate(name = factor(name,levels = c("TA","pH","Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL","Humics","Proteinaceous","HIX","Lignin_Like","M_C","MarineHumic_Like","VisibleHumidic_Like","Tryptophan_Like","Tyrosine_Like","TempInSitu_seep"))) %>%
  filter(!name %in% c("HIX","Lignin_Like","M_C","MarineHumic_Like","VisibleHumidic_Like","Tryptophan_Like","Tyrosine_Like") )%>%
  ggplot(aes(x = fct_reorder(nicenames, as.numeric(name), .desc = TRUE),#y = 1,
             y = c(rep(0.5,18),rep(1.5,18)),
             #y = c(rep(1,24),rep(1.8,24)),
             fill = -estimate, size = abs(estimate)))+
  geom_point(shape = 21)+
  scale_fill_gradient2(low = "#AD9024", high = "#AC7299",mid = "beige", midpoint = 0, limits = c(-1,1))+
  scale_size(range = c(0.1,10))+
  coord_flip()+
  theme_bw()+
  guides(size = "none",
         shape = "none",
         fill = guide_colourbar(title.position="top", title.hjust = 0.5))+
  labs(x = "",
       y = "",
       fill = "Correlation with freshwater"
         )+
  scale_y_continuous(breaks = c(.5,1.5), labels = c("Dry Season","Wet Season"), limits = c(0,2))+
                     #, limits = c(0.8,2.5))+
  theme( panel.grid = element_blank(),
        axis.text.y = element_markdown(size = 14),
        axis.text.x = element_markdown(size = 14),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 12),
         legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
      strip.background = element_blank(),
      strip.text = element_markdown(size = 14, face = "bold")
        )+
  facet_wrap(~Location, ncol = 2)

ggsave(here("Output","CorrelationPlot_seepSalinity.pdf"), height = 8, width = 8)



### Pure Spring Water
mean_spring<-Data %>%
  mutate(Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tyrosine_Like+Tryptophan_Like) %>%
  filter(CowTagID %in% c("VSPRING","Varari_Well","CSPRING_BEACH2","CSPRING_ROAD")) %>% # just do the seep data
  select(Location, Salinity,CowTagID, TA,pH,Phosphate_umolL:Ammonia_umolL, Humics,Proteinaceous, TempInSitu_seep)%>%
  # select(Location, Salinity, TA: Lignin_Like, TempInSitu_seep) %>%
  pivot_longer(cols = c(Salinity, TA: TempInSitu_seep)) %>%
  drop_na()%>%
  group_by(CowTagID, name)%>%
  summarise(mean = round(min(value),2)) %>%
  mutate(unit = case_when(name == "TempInSitu_seep" ~ " &deg;C", # add units
                          name %in% c("Ammonia_umolL","NN_umolL","Silicate_umolL","Phosphate_umolL")~" &mu;mol L<sup>-1</sup>",
                          name %in% c("M_C","HIX")~ " Ratio",
                          name %in% c("Salinity")~ " psu",
                          name %in% c("pH")~ " pH total scale",
                          name == "TA"~ " &mu;mol kg<sup>-1</sup>",
                          name %in% c("Lignin_Like","Tyrosine_Like","Tryptophan_Like","MarineHumic_Like","VisibleHumidic_Like", "Humics","Proteinaceous")~" Raman Units"),
         nicenames = case_when(name == "TempInSitu_seep" ~ "Temperature",
                               name == "pH" ~ "pH<sub>T</sub>",
                               #    name == "Lignin_Like" ~"Lignin Like",
                               #   name == "M_C" ~ "M:C",
                               name == "Tyrosine_Like" ~ "Tyrosine Like",
                               name == "Tryptophan_Like" ~ "Tryptophan Like",
                               name == "HIX"~"HIX",
                               #  name == "MarineHumic_Like" ~ "Marine Humic Like",
                               name == "VisibleHumidic_Like" ~ "Visible Humic Like",
                               name == "Ammonia_umolL" ~ "Ammonium",
                               name == "TA" ~ "Total Alkalinity",
                               name == "Phosphate_umolL" ~ "Phosphate",
                               name == "NN_umolL" ~ "Nitrate+Nitrite",
                               name == "Silicate_umolL" ~ "Silicate",
                               name == "Salinity" ~"Salinity",
                               name == "Humics" ~"Humics fDOM",
                               name == "Proteinaceous"~"Proteinaceous fDOM"
                               
                               
         )) %>%
  pivot_wider(values_from = mean, names_from  = CowTagID)



mx<-mean_spring %>%
  mutate(name = factor(name,levels = c("Salinity","TA","pH","Ammonia_umolL","NN_umolL","Phosphate_umolL","Silicate_umolL","Humics","Proteinaceous","HIX","Lignin_Like","M_C","MarineHumic_Like","VisibleHumidic_Like","Tryptophan_Like","Tyrosine_Like","TempInSitu_seep"))) %>%
  column_to_rownames(var = "nicenames") %>%
  arrange(name)%>%
  select(unit,CSPRING_BEACH2,CSPRING_ROAD,VSPRING,Varari_Well ) %>%
  rename("Cabral Beach"=CSPRING_BEACH2, "Cabral Road"=CSPRING_ROAD,
         "Varari Beach"=VSPRING, "Varari Well"=Varari_Well) %>%
    htmlTable(css.cell = "width: 110px")

mx
   

### Only keep certain dataframes for eco metab script
rm(list= ls()[!(ls() %in% c("Data", "Datalog", "remove_varari","remove_cabral","remove_vararilog","remove_cabrallog"))])
