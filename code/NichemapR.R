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


options(timeout = 1000000000) #Stop nichemapr from timing out

#Run NicheMapR

#### 1: Load microclimate data   
df <- read.csv("data/NichmapR/NicheMapR dataset (Floristics).csv")

#### 2. Set dates
dstart <- "01/11/2021" #set the start date
dfinish <- "01/03/2022" #set the end date

### 3. Input to model
#Run loop through all columns
plan(multisession, workers = 4)
result_list <- future.apply::future_lapply(future.seed = NULL, 1:nrow(df), function(x) {
  nm <- NicheMapR::micro_ncep(loc = c(df$Long[x], df$Lat[x]), 
                              dstart = dstart, 
                              dfinish = dfinish, 
                              run.gads = 2,
                              solonly =1, #run if you only want SR
                              runshade = 0, #Assumed 0% shade on summits
                              Usrhyt = 0.5, #Local height (m) at which air temperature, wind speed and humidity are to be computed 
                              slope = df$Slope[x], aspect = df$Degree[x], 
                              IR = 0, #Clear-sky longwave radiation computed using Campbell and Norman(1998) 
                              maxshade = 5) 
  out <- as.data.frame(nm$metout)
  # Add site name to data
  out$Site <- df$Site[x]
  out$Aspect <- df$Aspect[x]
  out
})  %>%  dplyr::bind_rows()

##### calculate average daily SR for site between 8 and 6 pm for summer
#TIME = time of day (mins) in UTC time
#60min is UTC time -> to convert need to + 11 hrs (660min)
unique(result_list$TIME)

nichemap2 <- result_list %>% 
  filter(TIME <= "300") #filter all time before 4am (300min)
nichemap3 <- nichemap2 %>% 
  filter(TIME >= "1260")#then filter by times after 8am (1260min)
nichemap4 <- nichemap3[,c(1,2,13,20,21)] #select date, time, SR, site, aspect


#Calcuate the cumulative sum of SR for each aspect and each site 
#during the growing season

nichemap5 <- nichemap4 %>% 
  group_by(Site, Aspect) %>% 
  summarise(avg.SR = sum(SOLR)) 


#Save as csv
#write.csv(nichemap5, "data/Cumulative SR for site per aspect.csv", row.names = FALSE)


###### END #######


