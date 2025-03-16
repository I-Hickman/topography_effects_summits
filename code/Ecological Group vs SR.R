
####### ECOLOGICAL GROUPS #####
### Endemic, generalist, exotic species cover & richness
#### Q: how does the cover of endemic, generalist and exotic species change with SR?

library(ggplot2)
library(cowplot)
library(brms)
library(dplyr)
library(gridExtra)
library(tibble)
library(tidyverse)
library(bayestestR)
library(GGally)
library(future)


## Step 1) Load data ######

Florisitc <- read.csv("data/2022 Summit florisitc.csv")
Traits <- read.csv("data/Trait dataset.csv")
SR <- read.csv("data/Cumulative SR for site per aspect (floristics).csv")

#Join SR to florisitics
Flor_SR <- left_join(Florisitc, SR, by = c("Site", "Aspect"))

#Make Florisitc long form
Florisitc_LO <- Flor_SR %>% 
  group_by(Year, Site, avg.SR, Quadrat) %>% 
  gather(key = "Species", value = "Cover", 14:122) 

#Join data sets
Flor_traits <- left_join(Florisitc_LO, Traits, by = "Species")

#remove any NAs in data (only Agrostis sp.)
Flor_traits <- Flor_traits %>% filter(!is.na(Origin))

#Check for any NA's in the origin column
anyNA(Flor_traits$Origin)

#Remove SR.w.m2 column
Flor_traits <- Flor_traits %>% select(-SR.w.m2)

#Examine data for any issues
head(Flor_traits)
unique(Flor_traits$avg.SR)

## Step 2) Explore data ######
#Find the average frequency of species occurrence per quadrat
average_frequency <- Florisitc_LO %>%
  filter(Cover != 0) %>%
  group_by(Year, Species) %>%
  summarise(total_occurrences = n()) 
  
average_frequency2 <- average_frequency %>% 
  group_by(Year) %>%
  summarise(frequency = mean(total_occurrences))

#highest frequency species
#Carex brevi = 205
#number of quadrats = 280
205/280

#Find out the number of endemic, native and exotic species
count<- Flor_traits %>%
  group_by(Species, Growth_form) %>%
  summarise(n = n())

count %>%
  group_by(Growth_form) %>%
  summarise(n = n())

#Higher cover
Flor_traits %>%
  group_by(Species) %>%
  summarise(total_cover = sum(Cover)) %>%
  arrange(-total_cover) %>%
  head(10)
#What total percent of cover did Poa.fawcettiae contribute to the total cover?
sum(Flor_traits$Cover[Flor_traits$Species == "Poa.fawcettiae"])/sum(Flor_traits$Cover)
#What total percent of cover did Hovea.montana  contribute to the total cover?
sum(Flor_traits$Cover[Flor_traits$Species == "Hovea.montana"])/sum(Flor_traits$Cover)




############ 1) RICHNESS  ###########

## Step 1. Calculate richness for each quadrat ######
Flor_richness <-Flor_traits %>%
  group_by(Year, Site, Aspect, ele, avg.SR, Quadrat, Origin) %>%
  filter(Cover != "0") %>% 
  summarise(species.richness=n()) %>%
  arrange(-species.richness)


## Step 2) Model #####
## Scale the independent variable (SR)
Flor_richness$SR.std = as.vector(scale(Flor_richness$avg.SR)) 

## Create model
mRich <- brms::brm(species.richness ~ SR.std*Origin + (1 | Site),
                   family = poisson,
                   iter = 4000,
                   cores = 4,
                   chains = 4,
                   seed = 1234,
                   control = list(adapt_delta = 0.8, max_treedepth=15),
                   data = Flor_richness)

## Step 3) Validate model #######
#Check model
summary(mRich)
plot(mRich)
pp_check(mRich) 


## Step 4) Fixed effects coefficient plot #####
mRich_fixed <- as.data.frame(fixef(mRich))
mRich_fixed2 <- rownames_to_column(mRich_fixed, var = "Variable")
mRich_fixed2$Variable <- factor(mRich_fixed2$Variable,
                               levels = c("SR.std:OriginNative", 
                                          "SR.std:OriginExotic", 
                                          "OriginNative", 
                                          "OriginExotic", 
                                          "SR.std", 
                                          "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "SR.std:OriginNative" = "SR:Generalist",
                   "SR.std:OriginExotic" = "SR:Exotic",
                   "OriginNative" = "Generalist",
                   "OriginExotic" = "Exotic",
                   "Intercept" = "Intercept")

#Get the 80% credible intervals from model
ci_eti_fixed_95 <- bayestestR::ci(mRich, method = "ETI", ci = 0.95, effects = "fixed") #95% CI

#Rename CI low and CI high columns 
ci_eti_fixed_95 <- ci_eti_fixed_95 %>%
  rename(Q2.5 = CI_low, Q97.5 = CI_high) #95% CI

# add in estimate column 
ci_eti_fixed_95$Estimate <- mRich_fixed2$Estimate
ci_eti_fixed_95$Variable <- mRich_fixed2$Variable

End_fixed2 <- ggplot(ci_eti_fixed_95, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  #geom_errorbarh(aes(xmin = Q5.5, xmax = Q94.5, height = 0), linewidth = 0.7) +  # 89% Confidence intervals)
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0), linewidth = 0.3) +  # 95% Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")


## Step 5) Random effects ####

mRich_rand <- as.data.frame(ranef(mRich))
mRich_rand2 <- rownames_to_column(mRich_rand, var = "Variable")
mRich_rand2 <- mRich_rand2[order(mRich_rand2$Site.Estimate.Intercept), ]

End_rand <- ggplot(mRich_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Endemic, gen & exotics/Richness/EGE richness model (random effects).png", 
       plot = End_rand, width = 10, height = 7.07) 


## Step 6) Plot predictions ####

mean(Flor_richness$avg.SR) #461761.6
sd(Flor_richness$avg.SR) #20812.24
range(Flor_richness$avg.SR) #393523.9 486082.5

#ENDEMIC SPECIES
new_end1 <- data.frame(Site = "average site", 
                      Origin = "Endemic",
                      SR.w.m2 = seq(393520, 486090, by = 10)) %>% 
  dplyr::mutate(SR.std = (SR.w.m2 - 461761.6)/20812.24)
#Predictions
predictions<- posterior_epred(mRich, newdata = new_end1,  
                              allow_new_levels = TRUE) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_end1$mean <- mean
new_end1$lci <- lci
new_end1$uci <- uci

#GENERALIST SPECIES
new_gen1 <- data.frame(Site = "average site", 
                      Origin = "Native",
                      SR.w.m2 = seq(393520, 486090, by = 10)) %>% 
  dplyr::mutate(SR.std = (SR.w.m2 - 461761.6)/20812.24)
#Predictions
predictions<- posterior_epred(mRich, newdata = new_gen1,  
                              allow_new_levels = TRUE) 

#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_gen1$mean <- mean
new_gen1$lci <- lci
new_gen1$uci <- uci

#EXOTIC SPECIES
new_exotic1 <- data.frame(Site = "average site", 
                      Origin = "Exotic",
                      SR.w.m2 = seq(393520, 486090, by = 10)) %>% 
  dplyr::mutate(SR.std = (SR.w.m2 - 461761.6)/20812.24)
#Predictions
predictions<- posterior_epred(mRich, newdata = new_exotic1,  
                              allow_new_levels = TRUE) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_exotic1$mean <- mean
new_exotic1$lci <- lci
new_exotic1$uci <- uci


### Plot
EGE_plot <- ggplot() + 
  # Endemic species
  geom_path(data = new_end1, aes(x = SR.w.m2, y = mean, color = "Endemic")) +
  geom_ribbon(data = new_end1, aes(x = SR.w.m2, ymin = lci, ymax = uci), alpha = 0.2, fill = "darkgreen") +
  # Generalist species
  geom_path(data = new_gen1, aes(x = SR.w.m2, y = mean, color = "Native")) +
  geom_ribbon(data = new_gen1, aes(x = SR.w.m2, ymin = lci, ymax = uci), alpha = 0.2, fill = "darkblue") +
  # Exotic species
  geom_path(data = new_exotic1, aes(x = SR.w.m2, y = mean, color = "Exotic")) +
  geom_ribbon(data = new_exotic1, aes(x = SR.w.m2, ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  theme_cowplot() +
  ylab(expression("Species richness")) +
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) +
  scale_color_manual(values = c("Endemic" = "darkgreen", "Native" = "darkblue", "Exotic" = "darkorange")) +
  theme(legend.position = "right") +
  labs(color = "Ecological group")


#Add labels
End_fixed2 <- End_fixed2+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
EGE_plot <- EGE_plot+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)


## Step 7) Plot ####
plot_end <- grid.arrange(End_fixed2, EGE_plot, ncol = 2)

ggsave("output/Endemic, gen & exotics/Richness/EGE richness vs solar radiation.png", 
       plot = plot_end, width = 10, height = 3.5)




############## 2) PRESENCE & ABSENCES ############

#Explore data
hist(Flor_traits$Cover)

## Step 1: Convert cover data to Presence and absence data #####
EG_PA_data <- Flor_traits %>% 
  dplyr::select(Code, Site, SR = avg.SR, ele, Aspect, Quadrat,
                Species, Cover, Origin) %>%
  dplyr::mutate(PA = ifelse(Cover > 0, 1, 0))
  
#Look at histogram of data
hist(EG_PA_data$PA) 
unique(EG_PA_data$PA)

### Scale predictors
EG_PA_data$SR.std = as.vector(scale(EG_PA_data$SR)) 

## Step 2) Models ######
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE, threads_per_chain = 1)

#PA model
EG_model <- brms::brm(PA ~ SR.std*Origin + (1 | Site) + (1|Species),
                      family = bernoulli(),
                      iter = 4000, #increase iterations if issues 
                      cores = 4,
                      #future = TRUE,
                      chains = 4,
                      seed = 1234,
                      control = list(adapt_delta = 0.8, max_treedepth=16),
                      data = EG_PA_data)

## Step 3) Validate model #####
#Check model
summary(EG_model)
plot(EG_model)
pp_check(EG_model) 

## Step 4) Fixed effects coefficient plot ####
mEG_cov_fixed <- as.data.frame(fixef(EG_model))
mEG_cov_fixed <- rownames_to_column(mEG_cov_fixed, var = "Variable")
mEG_cov_fixed$Variable <- factor(mEG_cov_fixed$Variable,
                                 levels = c("SR.std:OriginNative", 
                                            "SR.std:OriginExotic", 
                                            "OriginNative", 
                                            "OriginExotic", 
                                            "SR.std", 
                                            "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "SR.std:OriginNative" = "SR:Generalist",
                   "SR.std:OriginExotic" = "SR:Exotic",
                   "OriginNative" = "Generalist",
                   "OriginExotic" = "Exotic",
                   "Intercept" = "Intercept")

#Get the 95% credible intervals from model
ci_eti_fixed_95 <- bayestestR::ci(EG_model, method = "ETI", ci = 0.95, effects = "fixed") #95% CI

#Rename CI low and CI high columns 
ci_eti_fixed_95 <- ci_eti_fixed_95 %>%
  rename(Q2.5 = CI_low, Q97.5 = CI_high) #95% CI

#Add in estimate column 
ci_eti_fixed_95$Estimate <- mEG_cov_fixed$Estimate
ci_eti_fixed_95$Variable <- mEG_cov_fixed$Variable

#Plot
End_fixed <- ggplot(ci_eti_fixed_95, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  #geom_errorbarh(aes(xmin = Q5.5, xmax = Q94.5, height = 0), linewidth = 0.7) +  # 89% Confidence intervals)
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0), linewidth = 0.3) +  # 95% Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")



## Step 5) Random effects #####

# 1. Plot estimates of the standard deviation and CI of random effects to see which one explains the most in the unexplained variation 
summary(EG_model)
ran_DF <- data.frame(
  variable = c("Site", "Species"),
  estimate = c(0.37, 1.72),
  lci = c(0.25 ,0.58),
  uci = c(1.47, 2.00)
)

End_random_eff <- ggplot(ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")


### Site random effects
#SITE
mEG_rand <- as.data.frame(ranef(EG_model, groups = "Site"))
mEG_rand <- rownames_to_column(mEG_rand, var = "Variable")
mEG_rand <- mEG_rand[order(mEG_rand$Site.Estimate.Intercept), ]

EG_rand_St <- ggplot(mEG_rand, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

#SPECIES
mEG_rand <- as.data.frame(ranef(EG_model, groups = "Species"))
mEG_rand <- rownames_to_column(mEG_rand, var = "Variable")
mEG_rand <- mEG_rand[order(mEG_rand$Species.Estimate.Intercept), ]

EG_rand_Sp <- ggplot(mEG_rand, aes(x = Species.Estimate.Intercept, y = reorder(Variable, Species.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Species.Q2.5.Intercept, xmax = Species.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")



ggsave("output/Endemic, gen & exotics/Composition/EGE PA (Site random effects).png", 
       plot = EG_rand_St,width = 10, height = 7.07) 
ggsave("output/Endemic, gen & exotics/Composition/EGE PA (Species random effects).png", 
       plot = EG_rand_Sp, width = 10, height = 7.07) 


## Step 6) Plot predictions #########

mean(EG_PA_data$SR) #462020.3
sd(EG_PA_data$SR) #20579.3
range(EG_PA_data$SR) #393523.9 486082.5

#ENDEMIC SPECIES
new_end <- data.frame(Site = "average site", 
                      Species = "average species", 
                      Origin = "Endemic",
                      SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- posterior_epred(EG_model, newdata = new_end,  
                                     allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_end$mean <- mean
new_end$lci <- lci
new_end$uci <- uci

#GENERALIST SPECIES
new_gen <- data.frame(Site = "average site", 
                      Species = "average species", 
                      Origin = "Native",
                      SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- posterior_epred(EG_model, newdata = new_gen,  
                                     allow_new_levels = TRUE)

#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_gen$mean <- mean
new_gen$lci <- lci
new_gen$uci <- uci

#EXOTIC SPECIES
new_exotic <- data.frame(Site = "average site", 
                         Species = "average species", 
                         Origin = "Exotic",
                         SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- posterior_epred(EG_model, newdata = new_exotic,  
                                     allow_new_levels = TRUE)
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_exotic$mean <- mean
new_exotic$lci <- lci
new_exotic$uci <- uci


### Plot Endemics
End_cov_plot <- ggplot(new_end, aes(x = SR, y = mean)) + 
  geom_path(color = "darkgreen") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkgreen") +
  ylim(0,1) +
  theme_cowplot()+
  ylab(expression("Probability of occurance"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

### Plot Generalist
Gen_cov_plot <- ggplot(new_gen, aes(x = SR, y = mean)) + 
  geom_path(color = "darkblue") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkblue") +
  ylim(0,1) +
  theme_cowplot()+
  ylab(expression(" "))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

### Plot Exotics
Exot_cov_plot <- ggplot(new_exotic, aes(x = SR, y = mean)) + 
  geom_path(color = "darkorange") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  ylim(0,1) +
  theme_cowplot()+
  ylab(expression("Probability of occurance"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

#Add labels
End_fixed <- End_fixed+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
End_random_eff <- End_random_eff+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)
End_cov_plot <- End_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
Gen_cov_plot <- Gen_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)
Exot_cov_plot <- Exot_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)

## Step 7) Plot ######
plot_end <- grid.arrange(End_fixed, 
                         End_random_eff,
                         End_cov_plot,
                         Gen_cov_plot,
                         Exot_cov_plot, ncol = 2)

ggsave("output/Endemic, gen & exotics/Composition/EGE PA vs solar radiation3.png", 
       plot = plot_end, width = 11, height = 12)






############## 3) COVER ############

#Explore data
hist(Flor_traits$Cover)

## Step 1: Transform data #####
##Remove zeros
EG_cov_data <- Flor_traits %>% 
  dplyr::filter(Cover > 0) %>%
  dplyr::select(Code, Site, SR = avg.SR, ele, Aspect, Quadrat,
                Species, Cover, Origin) %>%
  dplyr::mutate(pCover = Cover/100) #turn % cover into proportional data
EG_cov_data <- EG_cov_data %>% 
  dplyr::mutate(qCover = qlogis(pCover))

anyNA(EG_cov_data)

#Look at histogram of data
hist(EG_cov_data$pCover)
hist(EG_cov_data$qCover)
unique(EG_cov_data$qCover)

### Scale predictors
EG_cov_data$SR.std = as.vector(scale(EG_cov_data$SR)) 

## Step 2) Models ######
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE, threads_per_chain = 1)

#Cover model
EG_cov_model <- brms::brm(qCover ~ SR.std*Origin + (1 | Site) + (1|Species),
                         family = gaussian(),
                         iter = 4000, #increase iterations if issues 
                         cores = 4,
                         #future = TRUE,
                         chains = 4,
                         control = list(adapt_delta = 0.8, max_treedepth=16),
                         data = EG_cov_data)

## Step 3) Validate model #####
#Check model
summary(EG_cov_model)
plot(EG_cov_model)
pp_check(EG_cov_model) 

## Step 4) Fixed effects coefficient plot ####
mEG_cov_fixed <- as.data.frame(fixef(EG_cov_model))
mEG_cov_fixed <- rownames_to_column(mEG_cov_fixed, var = "Variable")
mEG_cov_fixed$Variable <- factor(mEG_cov_fixed$Variable,
                                 levels = c("SR.std:OriginNative", 
                                            "SR.std:OriginExotic", 
                                            "OriginNative", 
                                            "OriginExotic", 
                                            "SR.std", 
                                            "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "SR.std:OriginNative" = "SR:Generalist",
                   "SR.std:OriginExotic" = "SR:Exotic",
                   "OriginNative" = "Generalist",
                   "OriginExotic" = "Exotic",
                   "Intercept" = "Intercept")

#Get the 95% credible intervals from model
ci_eti_fixed_95 <- bayestestR::ci(EG_cov_model, method = "ETI", ci = 0.95, effects = "fixed") #95% CI

#Rename CI low and CI high columns 
ci_eti_fixed_95 <- ci_eti_fixed_95 %>%
  rename(Q2.5 = CI_low, Q97.5 = CI_high) #95% CI

#Add in estimate column 
ci_eti_fixed_95$Estimate <- mEG_cov_fixed$Estimate
ci_eti_fixed_95$Variable <- mEG_cov_fixed$Variable

#Plot
End_fixed <- ggplot(ci_eti_fixed_95, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  #geom_errorbarh(aes(xmin = Q5.5, xmax = Q94.5, height = 0), linewidth = 0.7) +  # 89% Confidence intervals)
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0), linewidth = 0.3) +  # 95% Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")



## Step 5) Random effects #####

# 1. Plot estimates of the standard deviation and CI of random effects to see which one explains the most in the unexplained variation 
summary(EG_cov_model)
Forbs_ran_DF <- data.frame(
  variable = c("Site", "Species"),
  estimate = c(0.38, 0.94),
  lci = c(0.25 ,0.59),
  uci = c(0.57, 1.11)
)

End_random_eff <- ggplot(Forbs_ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")


### Site random effects
#SITE
mEG_rand <- as.data.frame(ranef(EG_cov_model, groups = "Site"))
mEG_rand <- rownames_to_column(mEG_rand, var = "Variable")
mEG_rand <- mEG_rand[order(mEG_rand$Site.Estimate.Intercept), ]

EG_rand_St <- ggplot(mEG_rand, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

#SPECIES
mEG_rand <- as.data.frame(ranef(EG_cov_model, groups = "Species"))
mEG_rand <- rownames_to_column(mEG_rand, var = "Variable")
mEG_rand <- mEG_rand[order(mEG_rand$Species.Estimate.Intercept), ]

EG_rand_Sp <- ggplot(mEG_rand, aes(x = Species.Estimate.Intercept, y = reorder(Variable, Species.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Species.Q2.5.Intercept, xmax = Species.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")



ggsave("output/Endemic, gen & exotics/Composition/EGE cover (Site random effects).png", 
       plot = EG_rand_St,width = 10, height = 7.07) 
ggsave("output/Endemic, gen & exotics/Composition/EGE cover (Species random effects).png", 
       plot = EG_rand_Sp, width = 10, height = 7.07) 


## Step 6) Plot predictions #########

mean(EG_cov_data$SR) #462020.3
sd(EG_cov_data$SR) #20579.3
range(EG_cov_data$SR) #393523.9 486082.5

#ENDEMIC SPECIES
new_end <- data.frame(Site = "average site", 
                      Species = "average species", 
                      Origin = "Endemic",
                      SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- plogis(posterior_epred(EG_cov_model, newdata = new_end,  
                              allow_new_levels = TRUE)) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_end$mean <- mean
new_end$lci <- lci
new_end$uci <- uci

#GENERALIST SPECIES
new_gen <- data.frame(Site = "average site", 
                      Species = "average species", 
                      Origin = "Native",
                      SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- plogis(posterior_epred(EG_cov_model, newdata = new_gen,  
                              allow_new_levels = TRUE)) 

#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_gen$mean <- mean
new_gen$lci <- lci
new_gen$uci <- uci

#EXOTIC SPECIES
new_exotic <- data.frame(Site = "average site", 
                         Species = "average species", 
                         Origin = "Exotic",
                         SR = seq(393520, 486083, by = 100)) %>% 
  dplyr::mutate(SR.std = (SR - 462020.3)/20579.3)
#Predictions
predictions<- plogis(posterior_epred(EG_cov_model, newdata = new_exotic,  
                              allow_new_levels = TRUE)) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_exotic$mean <- mean
new_exotic$lci <- lci
new_exotic$uci <- uci


### Plot Endemics
End_cov_plot <- ggplot(new_end, aes(x = SR, y = mean)) + 
  geom_path(color = "darkgreen") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkgreen") +
  ylim(0,0.5) +
  theme_cowplot()+
  ylab(expression("Proportional cover"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

### Plot Generalist
Gen_cov_plot <- ggplot(new_gen, aes(x = SR, y = mean)) + 
  geom_path(color = "darkblue") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkblue") +
  ylim(0,0.5) +
  theme_cowplot()+
  ylab(expression(" "))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

### Plot Exotics
Exot_cov_plot <- ggplot(new_exotic, aes(x = SR, y = mean)) + 
  geom_path(color = "darkorange") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "darkorange") +
  ylim(0,0.5) +
  theme_cowplot()+
  ylab(expression("Proportional cover"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

#Add labels
End_fixed <- End_fixed+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
End_random_eff <- End_random_eff+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)
End_cov_plot <- End_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
Gen_cov_plot <- Gen_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)
Exot_cov_plot <- Exot_cov_plot+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)

## Step 7) Plot ######
plot_end <- grid.arrange(End_fixed, 
                         End_random_eff,
                         End_cov_plot,
                         Gen_cov_plot,
                         Exot_cov_plot, ncol = 2)

ggsave("output/Endemic, gen & exotics/Composition/EGE abundance vs solar radiation4.png", 
       plot = plot_end, width = 11, height = 12)




##### APPENDICES ######


#Richness
pdf("output/Richness model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mRich))
plot(preds$Estimate ~ mRich$data$species.richness,  
     xlab = "Richness observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 

#PA
pdf("output/EG PA model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(EG_model))
plot(preds$Estimate ~ EG_model$data$qCover,  
     xlab = "PA observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 


#Cover
pdf("output/EG Cover model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(EG_cov_model))
plot(preds$Estimate ~ EG_cov_model$data$qCover,  
     xlab = "Cover observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 
