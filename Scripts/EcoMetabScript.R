### Ecosystem metabolism script ####

####################################

# load libraries #################
library(tidyverse)
library(ggfortify)
library(lubridate)
library(ggforce)
library(viridis)
library(patchwork)
library(here)
library(ggridges)
library(seacarb)
library(ggtext)
library(lme4)
library(lmerTest)
library(broom)
library(oce)
library(margins)



# read in the data from the pca script so that all data cleaning is the same
source(here("Scripts","PCA_Script.R"))


## load the distance to seep and depth data to calculate residence times
distance<-read_csv(here("Data","distance_depth.csv")) %>%
  mutate(Depth_m = adj_CT_depth_cm/100)%>%
  select(CowTagID, dist_to_crest,Depth_m)


# calculate average depth of the site
SiteDepth<-distance %>%
  left_join(Data)%>%
  group_by(Location) %>%
  summarise(meanDepth = mean(Depth_m, na.rm = TRUE)) %>%
  mutate(reeflength = ifelse(Location  == "Varari",189,100)) # or 272 for cabral if south to north


## load the current speed data
currentspeed<-read_csv(here("Data","currentspeed.csv")) %>%
  mutate(DateTime = ymd_hms(DateTime))


#### Visualize the TA data ####

# bring both seasons together and add in other important metrics for calculating metabolism
Data <- Data %>%
  left_join(distance) %>% # add in the distance and depth data
  left_join(currentspeed)%>%
  left_join(SiteDepth) %>%
  select(!c(Top_Plate_ID:Jamie_Plate_ID)) %>%
  mutate(residence_time_h = (1/(CurrentSpeed_m_s*60))*reeflength,
          rho = swRho(salinity = Salinity, temperature = Temperature, pressure = 0, eos="unesco")) # calculate seawater density in kg/m3


### Calculate all carbonate parameters #############
#remove rows with NAs
Cdata<-Data %>%
   mutate(Tide_Time = paste(Tide, Day_Night)) %>%
  drop_na(TA) # keep only complete cases
 

#calculate rest of carbonate params
CO2<-carb(flag=8, Cdata$pH, Cdata$TA/1000000, S=Cdata$Salinity, 
          T=Cdata$Temperature, Patm=1, P=0, 
              k1k2="w14", kf="dg", ks="d", pHscale="T", b="u74", gas="potential")

#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propogation
er<-errors(flag=8, Cdata$pH, Cdata$TA/1000000, 
           S=Cdata$Salinity, T=Cdata$Temperature, 
           Patm=1, P=0,Pt=Cdata$Phosphate_umolL/1000000,
           Sit=Cdata$Silicate_umolL/1000000,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
sd(er$DIC*1000000)/sqrt(nrow(er))

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
CO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-CO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

Cdata[,c("CO2","HCO3","CO3","DIC","OmegaArag","OmegaCalcite","pCO2","fCO2")]<-
  CO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]


### Pull out the endmembers ###

### Calculate endmembers from the mixing line (remove the dug well and the street samples at Cabral)

Cdata %>% filter(Plate_Seep == "Seep"| Plate_Seep == "Spring") %>% filter( TA <5000 & TA >1500)

## calculate the endmembers from the springs
Vend_spring<-Cdata %>% # Varari
  filter(Plate_Seep == "Spring", CowTagID == "VSPRING") %>%
  select(Location,  TA_Spring = TA,DIC_Spring = DIC, Salinity_Spring = Salinity, Silicate_umolL_Spring = Silicate_umolL)%>%
  slice(1)

#cabral
Cend_spring<-Cdata %>% 
  filter(Plate_Seep == "Spring", Location == "Cabral", CowTagID %in% c("CSPRING_ROAD") )%>% # was CSPRING_ROAD prior
  select(
    Location,
    TA_Spring = TA, DIC_Spring = DIC, Salinity_Spring = Salinity, Silicate_umolL_Spring = Silicate_umolL)%>%
  slice(1)

Spring<-bind_rows(Vend_spring, Cend_spring)
## Calculate endmembers for the "ocean" which will be the average salinity or silicate during high tide by season

endmembers<- Cdata %>%
  filter(Plate_Seep == "Plate")%>%
  drop_na(Silicate_umolL)%>%
   group_by(Location) %>%
   summarise_at(vars(Salinity, Silicate_umolL), .funs = list(min, max)) %>%
   rename(Salinity_High = Salinity_fn2, Silicate_umolL_High = Silicate_umolL_fn1) %>%
   select(Location, Salinity_High, Silicate_umolL_High) %>%
   left_join(Spring) %>% # bring in the spring data
  mutate(Salinity_Ocean = 35.972,
         Silicate_umolL_Ocean = 0.965)


Cdata <- left_join(Cdata, endmembers)  # add the endmemebers to the dataframe to make it easier to calculate
  
## put them in the dataframe
Cdata <- Cdata %>% # add the predicted mixing line
    mutate(TA.mix = TA+(TA - TA_Spring)*((Silicate_umolL - .965)/(Silicate_umolL_Spring - Silicate_umolL)), # account for SGD mixing
           DIC.mix = DIC+(DIC - DIC_Spring)*((Silicate_umolL - .965)/(Silicate_umolL_Spring - Silicate_umolL)) )

## calculate average TA and DIC from offshore water/ VRC and CRS
crestwater<-Cdata %>%
  filter(Plate_Seep == "Offshore", DIC<2060) %>%
  group_by(Location, TimeBlock) %>%
  summarise(TA.offshore = mean(TA, na.rm = TRUE),
            DIC.offshore = mean(DIC, na.rm = TRUE),
            Salinity.offshore = mean(Salinity, na.rm = TRUE),
            Silicate_umolL.offshore = mean(Silicate_umolL, na.rm = TRUE)) %>%
  drop_na()

# Data from Global Ocean Data Project
#GLODAPv2 Merged Dataset GLODAP Olsen, A., R. M. Key, S. van Heuven, S. K. Lauvset, A. Velo, X. Lin, C. Schirnick, A. Kozyr, T. Tanhua, M. Hoppema, S. Jutterström, R. Steinfeldt, E. Jeansson, M. Ishii, F. F. Pérez and T. Suzuki. The Global Ocean Data Analysis Project version 2 (GLODAPv2) – an internally consistent data product for the world ocean, Earth Syst. Sci. Data, 8, 297–323, 2016, (doi:10.5194/essd-8-297-2016).
carb(flag = 8, var1 = 8.134,var2 = 2364.613/1000000, S = 35.972, T = 25, Sit = .965/1000000, pHscale = "T", k1k2 = "x", b = "u74", kf = "dg", gas = "potential", ks = "d", P = 0)

OpenOcean<-tibble(TA_ocean = 2364.613, DIC_ocean = 1994.888, Si_ocean = 0.965, pH_ocean = 8.134, Temperature_ocean = 28.679)

# Calculate NEC and NEP proxy
Cdata<-Cdata %>%
  left_join(crestwater)%>%
        mutate(
          TA.diff = TA.offshore - TA.mix, #TA_incom --- take average from offshore water 2360
          DIC.diff = DIC.offshore - DIC.mix, # DIC_income 2000
          NEC.proxy = TA.diff/2,
          NEP.proxy = DIC.diff - (TA.diff/2)
         )

