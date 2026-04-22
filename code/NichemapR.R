#################################
###### NicheMapR ############
############################

##Load packages

library(usethis)
library(devtools)
library(microclima)
library(NicheMapR)
library(dplyr)
library(future)
library(maps)
library(RNCEP)
library(future.apply)
library(microclima)

#Run NicheMapR

#### 1: Load microclimate data   
df <- read.csv("data/NicheMapR dataset (Floristics).csv")

#### 2. Set dates
dstart <- "01/11/2021" #set the start date
dfinish <- "01/03/2022" #set the end date

### 3. Input to model
#Run loop through all columns
options(timeout = 1000000000) #Stop nichemapr from timing out
plan(multisession, workers = 4)
result_list <- future.apply::future_lapply(future.seed = NULL, 1:nrow(df), function(x) {
  nm <- NicheMapR::micro_ncep(loc = c(df$Long[x], df$Lat[x]), 
                              dstart = dstart, 
                              dfinish = dfinish, 
                              run.gads = 2,
                              solonly =1, #run if you only want SR
                              runshade = 0, #Assumed 0% shade on summits
                              Usrhyt = 0.1, #Local height (m) at which air temperature, wind speed and humidity are to be computed 
                              slope = df$Slope[x], aspect = df$Degree[x], 
                              IR = 0, #Clear-sky longwave radiation computed using Campbell and Norman(1998) 
                              maxshade = 5) 
  out <- as.data.frame(nm$metout)
  # Add site name to data
  out$Site <- df$Site[x]
  out$Aspect <- df$Aspect[x]
  out
})  %>%  dplyr::bind_rows()

#4: Filter for daylight hours (8am-6pm local time)
#TIME is in minutes (UTC)
#Adjust for local timezone (AEDT = UTC+11, so add 660 minutes)
#8am local = 480 + 660 = 1140 min
#6pm local = 1080 + 660 = 1740 min

##### Calculate average daily SR for site between 8am and 6pm for summer
result_list <- result_list %>%
  mutate(
    # TIME is minutes from midnight UTC
    # Convert to hours for easier interpretation
    hour_UTC = TIME / 60,
    # Add 11 hours for AEDT (wrapping at 24)
    hour_local = (hour_UTC + 11) %% 24
  )

##### Calculate average daily SR for site between 8am and 6pm for summer
nichemap_daylight <- result_list %>% 
  group_by(DOY, Site, Aspect) %>%
  filter(hour_local >= 8 & hour_local <= 18) %>%  # 8am to 6pm AEDT
  select(DOY, hour_local, SOLR, Site, Aspect)

# Calculate cumulative solar radiation
nichemap_summary <- nichemap_daylight %>% 
  group_by(Site, Aspect) %>% 
  summarise(cumulative_SR = sum(SOLR), .groups = "drop")

#Save as csv
write.csv(nichemap_summary, "data/Cumulative SR for site per aspect2.csv", row.names = FALSE)


###### END #######


