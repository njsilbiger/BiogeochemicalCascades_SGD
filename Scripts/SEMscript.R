### SEM script for Moorea reef data ###
### By Nyssa Silbiger ####
### Created on 12/29/2022 ####

#### Load libraries ############
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(ggrepel)
library(patchwork)
library(brms)
library(tidybayes)
library(bayesplot)
library(stringr)
library(broom)



#### Read in the data ####
# read in the data from the pca script so that all data cleaning is the same
source(here("Scripts","EcoMetabScript.R"))

#remove unnecessary dataframes for this analysis
rm(list= ls()[!(ls() %in% c("Cdata","Data", "Datalog", "turb_all"))])

Cdata <- Cdata%>%
  filter(Silicate_umolL <15, NEC.proxy<100) # remove an outlier

# extract the model data and log transform variables that need it
ModelData<-Cdata %>%
  filter(Plate_Seep == "Plate") %>%
  select(Location, Tide, Day.Night = Day_Night, TimeBlock, Season, pH, Salinity, Silicate_umolL, NEP.proxy, NEC.proxy, NN_umolL, Phosphate_umolL, Ammonia_umolL, Temperature, VisibleHumidic_Like, MarineHumic_Like, Tryptophan_Like, Tyrosine_Like, TA) %>%
  mutate(SilicateumolL = log(Silicate_umolL), # log transform all the nutrient data
         NNumolL = log(NN_umolL),
         PhosphateumolL = log(Phosphate_umolL),
         AmmoniaumolL = log(Ammonia_umolL),
         Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tryptophan_Like +Tyrosine_Like) %>% 
  mutate_at(vars("pH":"Proteinaceous"),  function(x) scale(x, center = TRUE)) 

### Create a plot of the raw data that went into the model
ModelData_notlog<-Cdata %>%
  filter(Plate_Seep == "Plate") %>%
  mutate(Humics = VisibleHumidic_Like+MarineHumic_Like,
         Proteinaceous = Tryptophan_Like +Tyrosine_Like)%>%
  select(Location, Tide, Day.Night = Day_Night, TimeBlock, Season, Salinity,Temperature, Silicate_umolL, NN_umolL, TA, pH, Humics,Proteinaceous, NEP.proxy, NEC.proxy) %>%
  pivot_longer(cols = c(Salinity:NEC.proxy), names_to = "Parameters", values_to = "Values") %>%
  mutate(Parameters = factor(Parameters, levels= c("Salinity","Temperature", "Silicate_umolL", "NN_umolL","TA", "pH", "Humics","Proteinaceous", "NEP.proxy", "NEC.proxy"))) %>%
  mutate(nicenames = case_when(
    Parameters == "Salinity" ~"a. Salinity <br> (PSU)",
    Parameters == "Temperature" ~ "b. Temperature <br> (\u00B0C)",
    Parameters == "pH" ~ "c. pH<sub>T</sub>",
    Parameters == "TA" ~ "d. Total Alkalinity <br> (&mu;mol kg<sup>-1</sup>)",
    Parameters == "NN_umolL" ~ "e. Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
    Parameters == "Silicate_umolL" ~ "f. Silicate <br> (&mu;mol L<sup>-1</sup>)",
    Parameters == "Humics" ~ "g. Humic fDOM <br> (RU)",
    Parameters == "Proteinaceous" ~"h. Proteinaceous fDOM <br> (RU)",
    Parameters == "NEP.proxy"~"i. NEP Potential <br> (&Delta; &mu;mol kg<sup>-1</sup>)",
    Parameters == "NEC.proxy"~"j. NEC Potential <br> (&Delta; &mu;mol kg<sup>-1</sup>)"))

Datameans<-ModelData_notlog %>%
  group_by(Location, Season, TimeBlock, Tide, nicenames)%>%
  summarise(mean_vals = mean(Values, na.rm = TRUE),
            se_vals = sd(Values, na.rm = TRUE)/sqrt(n()))

p_dry<-ModelData_notlog %>%
  filter(Season == "Dry")%>%
  ggplot(aes(color = Tide, shape = TimeBlock))+
  geom_jitter(aes(x = Location, y = Values, color = Tide, shape = TimeBlock), alpha = 0.2)+
  geom_point(data = Datameans %>% filter(Season =="Dry"), aes(x = Location, y = mean_vals), size = 2)+
  geom_errorbar(data = Datameans %>% filter(Season =="Dry"), aes(x = Location, ymin = mean_vals - se_vals, ymax = mean_vals+se_vals), width = 0.1)+
  scale_color_manual(values = c("#85556e","#91b793"))+
  labs(x = "", y = " Dry Season Values",
       shape = "Time Block")+
  facet_wrap(~nicenames, scales = "free", nrow = 2)+
  theme_bw()+
  theme(strip.text = ggtext::element_markdown(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

ModelData_notlog <- ModelData_notlog %>%
  mutate(nicenames = case_when(
    Parameters == "Salinity" ~"k. Salinity <br> (PSU)",
    Parameters == "Temperature" ~ "l. Temperature <br> (\u00B0C)",
    Parameters == "pH" ~ "m. pH<sub>T</sub>",
    Parameters == "TA" ~ "n. Total Alkalinity <br> (&mu;mol kg<sup>-1</sup>)",
    Parameters == "NN_umolL" ~ "o. Nitrate+Nitrite <br> (&mu;mol L<sup>-1</sup>)",
    Parameters == "Silicate_umolL" ~ "p. Silicate <br> (&mu;mol L<sup>-1</sup>)",
    Parameters == "Humics" ~ "q. Humic fDOM <br> (RU)",
    Parameters == "Proteinaceous" ~"r. Proteinaceous fDOM <br> (RU)",
    Parameters == "NEP.proxy"~"s. NEP Potential <br> (&Delta; &mu;mol kg<sup>-1</sup>)",
    Parameters == "NEC.proxy"~"t. NEC Potential <br> (&Delta; &mu;mol kg<sup>-1</sup>)"))

Datameans<-ModelData_notlog %>%
  group_by(Location, Season, TimeBlock, Tide, nicenames)%>%
  summarise(mean_vals = mean(Values, na.rm = TRUE),
            se_vals = sd(Values, na.rm = TRUE)/sqrt(n()))

p_wet<-ModelData_notlog %>%
  filter(Season == "Wet")%>%
  ggplot(aes(color = Tide, shape = TimeBlock))+
  geom_jitter(aes(x = Location, y = Values, color = Tide, shape = TimeBlock), alpha = 0.2)+
  geom_point(data = Datameans %>% filter(Season =="Wet"), aes(x = Location, y = mean_vals), size = 2)+
  geom_errorbar(data = Datameans %>% filter(Season =="Wet"), aes(x = Location, ymin = mean_vals - se_vals, ymax = mean_vals+se_vals), width = 0.1)+
  scale_color_manual(values = c("#85556e","#91b793"))+
  labs(x = "", y = " Wet Season Values",
       shape = "Time Block")+
  facet_wrap(~nicenames, scales = "free", nrow = 2)+
  theme_bw()+
  theme(strip.text = ggtext::element_markdown(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

p_dry/p_wet+plot_layout(guides = "collect")

ggsave(here("Output","Data_raw.pdf"), width = 12, height = 10, device = cairo_pdf)
### Create the Bayesian linear models ####

NNmod<-bf(NNumolL ~SilicateumolL, family = "student") ### N and P are highly correlated 
Pmod<-bf(PhosphateumolL ~SilicateumolL, family = "student")
NH4mod<-bf(AmmoniaumolL ~SilicateumolL, family = "student") # N and NH4 are weakly correlated
NEPmod<-bf(NEP.proxy~(Day.Night*NNumolL) +Temperature, family = "student") ## N or P or NH4?
pHmod<-bf(pH~ NEP.proxy*SilicateumolL, family = "student")
NECmod<-bf(NEC.proxy~pH+Temperature, family = "student")
Humicsmod<-bf(Humics~(Day.Night*NEP.proxy)+SilicateumolL, family = "student")
Protmod<-bf(Proteinaceous~NEC.proxy+SilicateumolL, family = "student")

# Function to run Bayesian SEM and make the posterior predictive checks and plot marginal effects
RunSEM<-function(site, season){
  # fit the brms model
  fit_brms <- brm(  NNmod+
                    NEPmod+
                    pHmod+
                    NECmod+
                    Humicsmod+
                    Protmod+
                    set_rescor(FALSE),
    data=ModelData[ModelData$Location == site & ModelData$Season == season,]
    ,cores=4, chains = 3)
  
  # calculate LOO (leave one out) diagnostics
  SGD_loo<-loo(fit_brms, reloo = TRUE) # looks good!
  
  # plot posterior predictive checks
    p2<-pp_check(fit_brms, resp="NNumolL") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("N+N")
  p3<-pp_check(fit_brms, resp="NEPproxy") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("NEP")

  p4<-pp_check(fit_brms, resp="pH") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("pH")
  
  p5<-pp_check(fit_brms, resp="NECproxy") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("NEC")
  
  p6<-pp_check(fit_brms, resp="Humics") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("Humics")
  
  p7<-pp_check(fit_brms, resp="Proteinaceous") +
    scale_color_manual(values=c("red", "black"))+
    ggtitle("Proteinaceous")

 
    p2+p3+p4+p5+p6+p7+
    plot_layout(guides = "collect") +
    plot_annotation(title = 'Posterior Predictive Checks', tag_levels = "A")
  ggsave(here("Output",paste(site,season, "Posteriorchecks.pdf")), width = 5, height = 5)
  
  # plot the results
#   
  ## Model 2
  # NN ~ Silicate
  R<-conditional_effects(fit_brms, "SilicateumolL", resp = "NNumolL",  resolution = 200)
  Rdata2<-R$NNumolL.NNumolL_SilicateumolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center"),
           lower = lower__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center"),
           upper = upper__*attr(ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    ) %>%
    mutate(modname = "NN~Silicate",
           Location = site,
           Season = season)
  
  R2<-Rdata2 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = exp(estimate)), lwd = 1)+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = NN_umolL)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    coord_trans(x = "log", y="log")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')
  
  
  # Model 3
  # NEP ~ DayNight*(Temperature + logNN)
  R<-conditional_effects(fit_brms, "NNumolL:Day.Night", resp = "NEPproxy", resolution = 200)
  Rdata3<-R$`NEPproxy.NEPproxy_NNumolL:Day.Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")) %>%
    mutate(NNumolL = NNumolL*attr (ModelData$NNumolL,"scaled:scale")+attr(ModelData$NNumolL,"scaled:center")
    ) %>% mutate(modname = "NEP~NN*Day/Night",
                 Location = site,
                 Season = season)
  
  R3<-Rdata3 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(NNumolL), y = estimate, lty = Day.Night),  lwd = 1)+
    geom_ribbon(aes(x = exp(NNumolL),ymin=lower, ymax=upper, lty = Day.Night), fill = "grey", alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NN_umolL, y = NEP.proxy, color = Temperature, shape = Day_Night)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
    coord_trans(x = "log")+
    scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none",
         shape = guide_legend(title = "Day/Night", title.position="top", title.hjust = 0.5)) 
 

   # NEP ~ DayNight*(Temperature + logNN)
  R<-conditional_effects(fit_brms, "Temperature", resp = "NEPproxy", resolution = 200)
  Rdata3a<-R$`NEPproxy.NEPproxy_Temperature` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")) %>%
    mutate(Temperature = Temperature*attr (ModelData$Temperature,"scaled:scale")+attr(ModelData$Temperature,"scaled:center")
    ) %>% mutate(modname = "NEP ~ Temperature",
                 Location = site,
                 Season = season)
  
  R3a<- Rdata3a %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = Temperature, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = Temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEP.proxy, color = NN_umolL)) +
    ylab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    scale_color_gradient(name = "N+N", trans = "log", breaks =c(0.1,.2, .3, 0.5,1), low = "lightgreen", high = "darkgreen")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
 # model 4
  # pH ~ NEP *log(Silicate)
  R<-conditional_effects(fit_brms, "NEP.proxy", resp = "pH", resolution = 200)
  Rdata4<-R$pH.pH_NEP.proxy %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           lower = lower__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           upper = upper__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")) %>%
    mutate(NEP.proxy = NEP.proxy*attr (ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")
    )%>% mutate(modname = "pH ~ NEP",
                Location = site,
                Season = season)
    
  R4<- Rdata4 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEP.proxy, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = NEP.proxy,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = pH, color = Silicate_umolL)) +
    scale_color_gradient(name = "Silicate", trans = "log", breaks =c(0.1,1,2,3,5,10,15), low = "lightblue", high = "darkblue")+
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab("pH")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
      guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # pH ~ NEP + log(Silicate)
  
R<-conditional_effects(fit_brms, "SilicateumolL:NEP.proxy", resp = "pH", resolution = 200, re_formula = NULL)

  Rdata4a<-R$`pH.pH_SilicateumolL:NEP.proxy` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           lower = lower__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center"),
           upper = upper__*attr(ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    ) %>%
    mutate(modname = "pH~Silicate*NEP",
           Location = site,
           Season = season)
  
  R4a<-Rdata4a%>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = estimate, lty = factor(NEP.proxy)), lwd = 1)+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=lower, ymax=upper, lty = factor(NEP.proxy)), fill = "grey", alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = pH, color = NEP.proxy)) +
    scale_color_gradient2(low = "#D8B365",
                          mid = "gray88",
                          high = "#097969",
                          midpoint = 0)+
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab("pH")+
    coord_trans(x="log")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none") 
  
  # Model 5
  # NEC ~ pH + Temperature

  R<-conditional_effects(fit_brms, "pH", resp = "NECproxy", resolution = 200)
  Rdata5<-R$NECproxy.NECproxy_pH %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")) %>%
    mutate(pH = pH*attr (ModelData$pH,"scaled:scale")+attr(ModelData$pH,"scaled:center")
    ) %>%
    mutate(modname = "NEC ~ pH",
           Location = site,
           Season = season)
    
  R5<-Rdata5 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = pH, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = pH,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = pH, y = NEC.proxy, color = Temperature)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab("pH")+
    scale_color_gradient(low = "blue", high = "red", name = "Temperature")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # NEC ~ pH + Temperature
  
  R<-conditional_effects(fit_brms, "Temperature", resp = "NECproxy", resolution = 200)
  Rdata5a<-R$NECproxy.NECproxy_Temperature %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           lower = lower__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center"),
           upper = upper__*attr(ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")) %>%
    mutate(Temperature = Temperature*attr (ModelData$Temperature,"scaled:scale")+attr(ModelData$Temperature,"scaled:center")
    ) %>%
    mutate(modname = "NEC ~ Temperature",
           Location = site,
           Season = season)
  
  R5a<-Rdata5a %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = Temperature, y = estimate), lwd = 1)+
    geom_ribbon(aes(x = Temperature,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Temperature, y = NEC.proxy, color = pH)) +
    ylab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    xlab(expression(atop("Temperature", paste("(",degree, "C)"))))+
    scale_color_gradient(low = "peachpuff", high = "sienna", name = expression("pH"[t]))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) 
  
  # Model 6
  # Humics~NEP.proxy*Day.Night+Silicate_umolL
  
  R<-conditional_effects(fit_brms, "NEP.proxy:Day.Night", resp = "Humics", resolution = 200)
  Rdata6<-R$`Humics.Humics_NEP.proxy:Day.Night` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(NEP.proxy = NEP.proxy*attr (ModelData$NEP.proxy,"scaled:scale")+attr(ModelData$NEP.proxy,"scaled:center")
    ) %>%
    mutate(modname = "Humics ~ NEP*Day/Night",
           Location = site,
           Season = season)
  
  R6<-Rdata6 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEP.proxy, y = estimate, lty = Day.Night), lwd = 1)+
    geom_ribbon(aes(x = NEP.proxy,ymin=lower, ymax=upper, lty = Day.Night),fill = "grey", alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEP.proxy, y = VisibleHumidic_Like +MarineHumic_Like, color = Silicate_umolL, shape = Day_Night)) +
    scale_color_gradient(name = "Silicate", trans = "log", breaks =c(0.1,1,2,3,5,10,15), low = "lightblue", high = "darkblue")+
    xlab(expression(atop("NEP", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none",
           shape = guide_legend(title = "Day/Night", title.position="top", title.hjust = 0.5)) 
  
  # Humics~NEP.proxy*Day.Night+Silicate_umolL
  R<-conditional_effects(fit_brms, "SilicateumolL", resp = "Humics", resolution = 200)
  Rdata6a<-R$`Humics.Humics_SilicateumolL` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           lower = lower__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center"),
           upper = upper__*attr(ModelData$Humics,"scaled:scale")+attr(ModelData$Humics,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    ) %>%
    mutate(modname = "Humics ~ Silicate",
           Location = site,
           Season = season)
  
  R6a<-Rdata6a %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = estimate), lwd = 1)+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=lower, ymax=upper),  alpha=0.3, fill = 'grey')+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = VisibleHumidic_Like +MarineHumic_Like, color = NEP.proxy)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Humics","(RU)")))+
    coord_trans(x="log")+
    scale_color_gradient2(low = "#D8B365",
                          mid = "gray88",
                          high = "#097969",
                          midpoint = 0)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
           lty = "none")  
 
  # Model 7
  # bf(Proteinaceous~NEC.proxy+Day.Night)
  R<-conditional_effects(fit_brms, "NEC.proxy", resp = "Proteinaceous", resolution = 200)
  
  Rdata7<-R$`Proteinaceous.Proteinaceous_NEC.proxy` %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           lower = lower__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           upper = upper__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center")) %>%
    mutate(NEC.proxy = NEC.proxy*attr (ModelData$NEC.proxy,"scaled:scale")+attr(ModelData$NEC.proxy,"scaled:center")
    ) %>% 
    mutate(modname = "Prot ~ NEC",
           Location = site,
           Season = season)
  
  R7<-Rdata7 %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = NEC.proxy, y = estimate),  lwd = 1)+
    geom_ribbon(aes(x = NEC.proxy,ymin=lower, ymax=upper),fill = "grey", linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = NEC.proxy, y = Tryptophan_Like +Tyrosine_Like, color = Silicate_umolL)) +
    scale_color_gradient(name = "Silicate", trans = "log", breaks =c(0.1,1,2,3,5,10,15), low = "lightblue", high = "darkblue")+
    xlab(expression(atop("NEC", paste("(",Delta," ", mu, "mol kg"^-1,")"))))+
    ylab(expression(atop("Proteinaceous","(RU)")))+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
  
  # Model 7
  R<-conditional_effects(fit_brms, "SilicateumolL", resp = "Proteinaceous", resolution = 200)
  Rdata7a<-R$Proteinaceous.Proteinaceous_SilicateumolL %>% # back transform the scaled effects for the plot
    mutate(estimate = estimate__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           lower = lower__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center"),
           upper = upper__*attr(ModelData$Proteinaceous,"scaled:scale")+attr(ModelData$Proteinaceous,"scaled:center")) %>%
    mutate(SilicateumolL = SilicateumolL*attr (ModelData$SilicateumolL,"scaled:scale")+attr(ModelData$SilicateumolL,"scaled:center")
    ) %>%
    mutate(modname ="Prot ~ Silicate",
           Location = site,
           Season = season)
  
  R7a<-Rdata7a %>%
    ggplot()+ # back trasform the log transformed data for better visual
    geom_line(aes(x = exp(SilicateumolL), y = estimate),  lwd = 1)+
    geom_ribbon(aes(x = exp(SilicateumolL),ymin=lower, ymax=upper),fill = "grey", linetype=1.5, alpha=0.3)+
    geom_point(data = Cdata %>% filter(Plate_Seep == "Plate", Location  == site, Season == season), aes(x = Silicate_umolL, y = Tryptophan_Like +Tyrosine_Like, color = NEC.proxy)) +
    xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
    ylab(expression(atop("Proteinaceous","(RU)")))+
    coord_trans(x="log")+
    scale_color_gradient2(low = "lightpink",
                          mid = "gray88",
                          high = "purple",
                          midpoint = 0)+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom',
          legend.box = "horizontal")+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
  
  R1<-ggplot() + # empty plot placeholder
    geom_text(aes(x = 1, y = 1, label = paste(site, season)), size = 10) +
    theme_void()+
    labs(x = "", 
         y = "")+
    scale_x_continuous(labels = NULL)+
    scale_y_continuous(labels = NULL)+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
    ## bring them all together in patchwork
  R<-(R1|R2)/(R3|R3a)/(R4|R4a)/(R5|R5a)/(R6|R6a)+(R7|R7a)+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = "collect")&
    guides(shape = guide_legend(order = 1, title.position = "top", title = "Day/Night",title.hjust = 0.5))&
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 11))
  
  ggsave(here("Output",paste(site,season,"marginaleffects.pdf")),R, width = 12, height = 18, useDingbats = FALSE)
  
  var_name <- paste(fit_brms, site,season, sep="_") # Construct the name
  var_name <- str_replace_all(var_name, " ", "")   # remove white space
 
  
  ### bring together the conditional effects data from the plots
  cond_data<-bind_rows(Rdata2, Rdata3, Rdata3a, Rdata4, Rdata4a, Rdata5, Rdata5a, Rdata6, Rdata6a, Rdata7, Rdata7a)
  
  return(list(fit_brms, cond_data))
}

# Run the SEM
set.seed(6) # set the seed so we get the same answer every time
V_Dry_fit<-RunSEM(site ="Varari", season = "Dry")
V_Wet_fit<-RunSEM(site ="Varari", season = "Wet")
C_Dry_fit<-RunSEM(site ="Cabral", season = "Dry")
C_Wet_fit<-RunSEM(site ="Cabral", season = "Wet")

## Get the posterior

get_post<-function(site,season, fit){
  
  post <- as_draws_df(fit[[1]])
  
  Cof<-post %>% 
    select(starts_with("b"),-ends_with("Intercept")) %>%
    gather() %>% 
    group_by(key)%>%
    median_hdci()%>%
    mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
    separate(col = key,into = c("b", "dependent", "independent"),sep = "_")%>% #loose the b and bring the values back together
    mutate(independent  = recode(independent, Day = "Day/Night", `NEP.proxy:Day` = "NEP.proxy:Day/Night")) %>%# mutate(key = paste(dependent, independent))%>%
    mutate(dependent = factor(dependent, levels = c("SilicateumolL","NNumolL","NEPproxy","pH","NECproxy","Humics","Proteinaceous")))%>%
    mutate(dependent = recode(dependent,"SilicateumolL" = "Silicate","NNumolL" = "NOx","NEPproxy" = "NEP","NECproxy" = "NEC")) %>%
    mutate(independent = factor(independent, levels = c("NEC.proxy","Day.NightNight:NEP.proxy:SilicateumolL","Day.NightNight:NEP.proxy","Day.NightNight:SilicateumolL","pH","NEP.proxy:SilicateumolL","NEP.proxy","Day.NightNight","Day.NightNight:NNumolL","NNumolL","Temperature","SilicateumolL","Salinity")))%>%
    mutate(independent = recode(independent, "NEC.proxy" = "NEC","Day.NightNight:NEP.proxy:SilicateumolL" = "Day/Night x NEP x Silicate","Day.NightNight:NEP.proxy" = "Day/Night x NEP","Day.NightNight:SilicateumolL" = "Day/Night x Silicate","NEP.proxy:SilicateumolL" = "NEP x Silicate","NEP.proxy" = "NEP","Day.NightNight" = "Day/Night","Day.NightNight:NNumolL" = "Day/Night x NOx","NNumolL" = "NOx","SilicateumolL"= "Silicate"))%>%
    mutate(Location = site, Season = season)
  
  return(Cof)
}

VDryCoF<-get_post(site = "Varari",season = "Dry", fit = V_Dry_fit)
VWerCoF<-get_post(site = "Varari",season = "Wet", fit = V_Wet_fit)
CDryCoF<-get_post(site = "Cabral",season = "Dry", fit = C_Dry_fit)
CWetCoF<-get_post(site = "Cabral",season = "Wet", fit = C_Wet_fit)

#calculate the coefficients
AllCoefs <- VDryCoF %>%
  bind_rows(VWerCoF,CDryCoF,CWetCoF)

write_csv(AllCoefs, here("Data","Coefficients.csv")) # export the coefficients


#Make the coefficients plot
CoefPlot<-AllCoefs%>%
  mutate(Season = recode(Season, Dry = "Dry Season", Wet = "Wet Season")) %>%
  ggplot(aes(x = value, y = independent, shape = sig, color = Location)) +  # note 
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point(size = 3)+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  scale_shape_manual(values = c(1,16))+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(
    x = NULL, y = NULL) +
  theme_bw() +
  guides(shape = FALSE)+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
  )+
  facet_grid(Season~dependent, scales = "free_y", space='free')
ggplot2::ggsave("Output/coefficientsAll.pdf", width = 12, height = 8, useDingbats = FALSE)


##### Make the marginal effects plots with just the lines for each model

Allfits<-bind_rows(V_Dry_fit[[2]], V_Wet_fit[[2]], C_Dry_fit[[2]], C_Wet_fit[[2]]) # extract the fits

NNSi<-Allfits %>%
  filter(modname=="NN~Silicate") %>%
  ggplot(aes(x = exp(SilicateumolL), y = exp(estimate), color = Location, lty = Season))+
  geom_line(size = 1.5)+
  coord_trans(x = "log", y = "log") +
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("Silicate", " (",mu,"mol L"^-1,")")),
       y = expression(paste("N+N", " (",mu,"mol L"^-1,")")))+
  theme_classic()

NEPNN<-Allfits %>%
  filter(modname=="NEP~NN*Day/Night") %>%
  ggplot(aes(x = exp(NNumolL), y = estimate, color = Location, lty = Season))+
  geom_hline(aes(yintercept = 0), color = "grey")+
  geom_line(size = 1.5)+
  coord_trans(x = "log") +
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(y = expression(paste("NEP", " (mmol kg"^-1,")")),
       x = expression(paste("N+N", " (",mu,"mol L"^-1,")")))+
  theme_classic()+
  facet_wrap(~Day.Night)


NEPTemp<-Allfits %>%
  filter(modname=="NEP ~ Temperature") %>%
  ggplot(aes(x = Temperature, y = estimate, color = Location, lty = Season))+
  geom_hline(aes(yintercept = 0), color = "grey")+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(y = expression(paste("NEP", " (mmol kg"^-1,")")),
       x = expression(paste("Temperature", " (",degree,"C)")))+
  theme_classic()

NEPpH<-Allfits %>%
  filter(modname=="pH ~ NEP") %>%
  ggplot(aes(x = NEP.proxy, y = estimate, color = Location, lty = Season))+
  #geom_hline(aes(yintercept = 0), color = "grey")+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("NEP", " (mmol kg"^-1,")")),
       y = "pH")+
  theme_classic()


SipH<-Allfits %>%
  filter(modname=="pH~Silicate*NEP") %>%
  mutate(effect2__ = as.numeric(as.character(effect2__)),
    effect2 = case_when(effect2__ > 0.69 ~ "High NEP",
                             effect2__ < 0.69 & effect2__ >= - 0.32 ~ "Low NEP",
                             effect2__ < -0.32 ~ "High Net Respiration"),
    effect2 = factor(effect2, levels = c("High NEP","Low NEP","High Net Respiration"))) %>%
  ggplot(aes(x = exp(SilicateumolL), y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("Silicate", " (",mu,"mol L"^-1,")")),
       y = "pH")+
  coord_trans(x = "log")+
  theme_classic()+
  facet_wrap(~effect2)

NECpH<-Allfits %>%
  filter(modname=="NEC ~ pH") %>%
  ggplot(aes(x = pH, y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(y = expression(paste("NEC", " (mmol kg"^-1,")")),
       x = "pH")+
  theme_classic()

NECTemp<-Allfits %>%
  filter(modname=="NEC ~ Temperature") %>%
  ggplot(aes(x = Temperature, y = estimate, color = Location, lty = Season))+
  geom_hline(aes(yintercept = 0), color = "grey")+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(y = expression(paste("NEC", " (mmol kg"^-1,")")),
       x = expression(paste("Temperature", " (",degree,"C)")))+
  theme_classic()

NEPHumics<-Allfits %>%
  filter(modname=="Humics ~ NEP*Day/Night") %>%
  ggplot(aes(x = NEP.proxy, y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("NEP", " (mmol kg"^-1,")")),
       y = "Humics (RU)")+
  theme_classic()+
  facet_wrap(~Day.Night)

SiHumics<-Allfits %>%
  filter(modname=="Humics ~ Silicate") %>%
  ggplot(aes(x = exp(SilicateumolL), y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  coord_trans(x = "log")+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("Silicate", " (",mu,"mol L"^-1,")")),
       y = "Humics (RU)")+
  theme_classic()

NECProt<-Allfits %>%
  filter(modname=="Prot ~ NEC") %>%
  ggplot(aes(x = NEC.proxy, y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  labs(x = expression(paste("NEC", " (mmol kg"^-1,")")),
       y = "Proteinaceaous (RU)")+
  theme_classic()

SiProt<-Allfits %>%
  filter(modname=="Prot ~ Silicate") %>%
  ggplot(aes(x = exp(SilicateumolL), y = estimate, color = Location, lty = Season))+
  geom_line(size = 1.5)+
  scale_color_manual(values = c("#16697A","#D64550"))+
  coord_trans(x = "log")+
  labs(x = expression(paste("Silicate", " (",mu,"mol L"^-1,")")),
       y = "Proteinaceaous (RU)")+
  theme_classic()

# make an empty spave for the composit plot
R1 <- ggplot()+
  theme_void()+
  labs(x = "", 
       y = "")+
  scale_x_continuous(labels = NULL)+
  scale_y_continuous(labels = NULL)+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# bring all the marginal effects lines plots together
AllPlot<-(CoefPlot|NNSi)/(NEPNN|NEPTemp)/(NEPpH|SipH)/(NECpH|NECTemp)/(NEPHumics|SiHumics)+(NECProt|SiProt)+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect")&
  labs(color = "",
       lty = "")&
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.background = element_rect(color = "black")
          )

ggsave(here("Output","marginaleffects_lines.pdf"), width = 12, height = 18, useDingbats = FALSE)


# model for silicate and salinity
summary(lm(data = Cdata %>%
  filter(Plate_Seep == c("Plate"), Location == "Varari"), log(Silicate_umolL)~Salinity*Season))

summary(lm(data = Cdata %>%
             filter(Plate_Seep == c("Plate"), Location == "Cabral"), log(Silicate_umolL)~Salinity*Season))

