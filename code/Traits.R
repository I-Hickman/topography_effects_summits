############################
#SUMMIT TRAITS and SOLAR RADIATION
############################


library(Rcpp)
library(brms)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(glmmTMB)


#####TABLE OF CONTENTS
#
##
###
#### 1) Load data and select variables
#### 2) SHRUBS - abundance, traits and SR
#### 3) GRAMINOIDS - abundance, traits and SR
#### 4) FORBS - abundance, traits and SR
###
##
#

####### 1) Load data and select variables ####
Trait_data <- read.csv("data/Trait and florisitic data.csv") #this has no 0's in it

#Omit species with missing trait data
Trait_data = Trait_data %>% na.omit()

#Plot the data and check for linear relationship
plot(qlogis(Trait_data$Cover/100) ~ log(Trait_data$avg.SR))
hist(Trait_data$Cover)
hist(Trait_data$avg.SR)

#Select variables & transform cover data
Trait_data2 <- Trait_data %>% 
  dplyr::select(Code, Site, SR =avg.SR, Quadrat,
                Species, Cover, Origin, Growth_form, 
                SLA = Mean_SLA_mm2.mg.1, LDMC = Mean_LDMC_mg.g.1, 
                Seed_mass = Mean_seed_mass_mg, Height = Max_height_m) %>%
  dplyr::mutate(pCover = (Cover/100)+0.001) #scaling (+0.001)

Trait_data2 <- Trait_data2 %>% 
  dplyr::mutate(qCover = qlogis(pCover)) #logit transformation

glimpse(Trait_data2) 

#Separate into lifeforms for analysis
Shrubs <- Trait_data2 %>% 
  filter(Growth_form == "Shrub")
Forbs <- Trait_data2 %>% 
  filter(Growth_form == "Forb")
Graminoid <- Trait_data2 %>% 
  filter(Growth_form == "Graminoid")




############## SHRUB SPECIES #################


glimpse(Shrubs)

## Step 1) Check correlation between variables #####
#Only put variables in the model that have biological meaning and this goes the same for interactions
Shrubs %>%
  dplyr::select(-c(Code, Site, Quadrat, Species, Origin, Growth_form, qCover, pCover, Cover)) %>%
  GGally::ggpairs()
#Low correlation between traits


## Step 2) Standardise fixed and random variables #####
#Scaled variables = (Shrubs$LDMC - mean)/sd(Shrubs$LDMC)
# I.e. Standardise variables (mean = 0 and SD = 1)
Shrubs$SR.std = as.vector(scale(Shrubs$SR))
Shrubs$SLA.std = as.vector(scale(Shrubs$SLA))
Shrubs$LDMC.std = as.vector(scale(Shrubs$LDMC))
Shrubs$Seed_mass.std = as.vector(scale(Shrubs$Seed_mass))
Shrubs$Height.std = as.vector(scale(Shrubs$Height))
glimpse(Shrubs)


## Step 3) Model #####
#Shrubs n=21
#Adding in variables/traits according the most importance
#most interested in SR, then LDMC, Height, then seed mass

mTraitsShrubs <- brms::brm(qCover ~
                    LDMC.std + SLA.std + Height.std + Seed_mass.std + SR.std +
                    LDMC.std:SR.std + Height.std:SR.std + Seed_mass.std:SR.std + 
                      SLA.std:SR.std +
                    (1|Site) + (1|Species), 
                    #could add (SR.std|)
                  family = gaussian,
                  iter = 4000,
                  cores = 4,
                  chains = 4,
                  seed= 1234,
                  control = list(adapt_delta = 0.9, max_treedepth=15),
                  data = Shrubs)

## Step 4) Validate model ####
#Check model
summary(mTraitsShrubs)
plot(mTraitsShrubs)
pp_check(mTraitsShrubs)
prior_summary(mTraitsShrubs)


## Step 5) Fixed effects coefficient plot #####
TrShrub_fixed <- as.data.frame(fixef(mTraitsShrubs))
TrShrub_fixed2 <- rownames_to_column(TrShrub_fixed, var = "Variable")
TrShrub_fixed2$Variable <- factor(TrShrub_fixed2$Variable,
                               levels = c("LDMC.std:SR.std","SLA.std:SR.std","Height.std:SR.std","Seed_mass.std:SR.std", 
                                          "LDMC.std", "SLA.std", "Height.std","Seed_mass.std","SR.std",
                                          "Intercept"))

custom_labels <- c("LDMC.std" = "LDMC", "SLA.std" = "SLA", "Height.std" = "Height", "Seed_mass.std" = "Seed mass", "SR.std" = "SR", 
                   "LDMC.std:SR.std" = "LDMC:SR", "SLA.std:SR.std" = "SLA:SR", "Height.std:SR.std" = "Height:SR", "Seed_mass.std:SR.std" = "Seed mass:SR",
                   "(Intercept)" = "Intercept")

TrShrub_fixed <- ggplot(TrShrub_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")

## Step 6) Random effects #####
# 1. Plot estimates of the standard deviation and CI of random effects to see which one explains the most in the unexplained variation 
summary(mTraitsShrubs)
shrub_ran_DF <- data.frame(
  variable = c("Site", "Species"),
  estimate = c(0.34, 0.76),
  lci = c(0.17, 0.43),
  uci = c(0.59, 1.22)
)

Shrub_rand <- ggplot(shrub_ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

ggsave("output/Traits/Bayesian/Shrubs/Shrub traits (random effect).png", 
       plot = Shrub_rand, width = 10, height = 7.07) 


# 2. Site
TrShrub_randSite <- as.data.frame(ranef(mTraitsShrubs, groups = "Site"))
TrShrub_randSite2 <- rownames_to_column(TrShrub_randSite, var = "Variable")
TrShrub_randSite2 <- TrShrub_randSite2[order(TrShrub_randSite2$Site.Estimate.Intercept), ]

Shrub_Site_rand <- ggplot(TrShrub_randSite2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")

ggsave("output/Traits/Bayesian/Shrubs/Shrub traits (Site random).png", 
       plot = Shrub_Site_rand, width = 10, height = 7.07) 

# 3. Species
TrShrub_randSpp <- as.data.frame(ranef(mTraitsShrubs, groups = "Species"))
TrShrub_randSpp2 <- rownames_to_column(TrShrub_randSpp, var = "Variable")
TrShrub_randSpp2 <- TrShrub_randSpp2[order(TrShrub_randSpp2$Species.Estimate.Intercept), ]

Shrub_Species_rand <- ggplot(TrShrub_randSpp2, aes(x = Species.Estimate.Intercept, y = reorder(Variable, Species.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Species.Q2.5.Intercept, xmax = Species.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Traits/Bayesian/Shrubs/Shrub traits (Species random).png", 
       plot = Shrub_Species_rand, width = 10, height = 7.07) 


## Step 7) Plot predictions ####
### LOW & HIGH SR
mean(Shrubs$SR) #463085.5
sd(Shrubs$SR) #19306.89
range(Shrubs$SR) #393523.9 486082.5
##LOW SR (at 393523.9W/m2 - back calcuated (393523.9 -mean)/SD = -2.285559)
(393523.9 - 463085.5)/19306.89 = -3.602942
#HIGH SR (at 486082.5W/m2 - back calcuated (486082.5 -mean)/SD = -2.285559)
(486082.5 - 463085.5)/19306.89 = 1.191129


##### LDMC #####
mean(Shrubs$LDMC) #491.035
sd(Shrubs$LDMC) #79.58581
range(Shrubs$LDMC) #296.1039 646.1538

#Create new dataset to predict to the average site when other variables at at 0 (their mean) and 
#### LOW LDMC ###
##1. Predict LDMC when SR is low 
newdat_shrub_ldmc_low <- data.frame(Site = "average site", #predicting to the average site
                              Species = "average species", #predicting to the average species
                              LDMC = seq(290, 650, by = 5), #range() of LDMC
                              Seed_mass.std = 0, Height.std = 0, SLA.std = 0, SR.std = -3.602942) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 491.035)/79.58581) #back scaling
newdat_shrub_ldmc_low_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrub_ldmc_low, 
                                                   allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_shrub_ldmc_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrub_ldmc_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrub_ldmc_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrub_ldmc_low$mean <- mean
newdat_shrub_ldmc_low$lci <- lci
newdat_shrub_ldmc_low$uci <- uci


###HIGH LDMC
#2. Predict LDMC when SR is high at (980 -918.8878)/38.89106 = 1.571369
newdat_shrub_ldmc_high <- data.frame(Site = "average site", #predicting to the average site
                               Species = "average species", #predicting to the average species
                               LDMC = seq(290, 650, by = 5), #range() of LDMC
                               Seed_mass.std = 0, Height.std = 0, SLA.std = 0, SR.std = 1.191129) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 491.035)/79.58581) #back scaling
newdat_shrub_ldmc_high_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrub_ldmc_high, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_shrub_ldmc_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrub_ldmc_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrub_ldmc_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrub_ldmc_high$mean <- mean
newdat_shrub_ldmc_high$lci <- lci
newdat_shrub_ldmc_high$uci <- uci



### LDMC PLOT
Shrub_LDMC_plot <- ggplot(newdat_shrub_ldmc_high, aes(x = LDMC, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_shrub_ldmc_low, aes(x = LDMC, y = mean), color = "blue") +
  geom_ribbon(data = newdat_shrub_ldmc_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression("Shrub species proportional cover")) +
  xlab(expression(LDMC ~ (mg/g^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))




##### Height ######
mean(Shrubs$Height) #0.7723404
sd(Shrubs$Height) #0.4732705
range(Shrubs$Height) #0.1 2.0

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_shrub_height_low <- data.frame(Site = "average site", #predicting to the average site
                                Species = "average species", #predicting to the average species
                                Height = seq(0.1, 2, by = 0.05), #range() of LDMC
                                LDMC.std = 0, Seed_mass.std = 0, SLA.std = 0, SR.std = -3.602942) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Height.std = (Height - 0.7723404)/0.4732705) #back scaling
newdat_height_low_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrub_height_low, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_height_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_height_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_height_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrub_height_low$mean <- mean
newdat_shrub_height_low$lci <- lci
newdat_shrub_height_low$uci <- uci



#HIGH SR
newdat_shrub_height_high <- data.frame(Site = "average site", #predicting to the average site
                                 Species = "average species", #predicting to the average species
                                 Height = seq(0.1, 2, by = 0.05), #range() of LDMC
                                 LDMC.std = 0, Seed_mass.std = 0, SLA.std = 0, SR.std = 1.191129) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Height.std = (Height - 0.7723404)/0.4732705) #back scaling)
newdat_height_high_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrub_height_high, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_height_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_height_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_height_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrub_height_high$mean <- mean
newdat_shrub_height_high$lci <- lci
newdat_shrub_height_high$uci <- uci



#Plot
Shrub_height_plot <- ggplot(newdat_shrub_height_high, aes(x = Height, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_shrub_height_low, aes(x = Height, y = mean), color = "blue") +
  geom_ribbon(data = newdat_shrub_height_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression("Shrub species proportional cover")) +
  xlab(expression("Height (m)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))



##### Seed mass #####
mean(Shrubs$Seed_mass) #8.32813
sd(Shrubs$Seed_mass) #6.19139
range(Shrubs$Seed_mass) # 0.0196 17.3861

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_shrubs_SM_low <- data.frame(Site = "average site", #predicting to the average site
                            Species = "average species", #predicting to the average species
                            Seed_mass = seq(0.0196, 18, by = 0.5), #range() of LDMC
                            LDMC.std = 0, Height.std = 0, SLA.std = 0, SR.std = -3.602942) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = (Seed_mass - 8.32813)/6.19139) #back scaling)
newdat_shrubs_SM_low_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrubs_SM_low, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_shrubs_SM_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrubs_SM_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrubs_SM_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrubs_SM_low$mean <- mean
newdat_shrubs_SM_low$lci <- lci
newdat_shrubs_SM_low$uci <- uci



##HIGH SR
newdat_shrubs_SM_high <- data.frame(Site = "average site", #predicting to the average site
                             Species = "average species", #predicting to the average species
                             Seed_mass = seq(0.0196, 18, by = 0.5), #range() of LDMC
                             LDMC.std = 0, Height.std = 0, SLA.std = 0, SR.std = 1.191129) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = (Seed_mass - 8.32813)/6.19139) #back scaling)
#Predictions
newdat_shrubs_SM_high_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrubs_SM_high, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_shrubs_SM_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrubs_SM_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrubs_SM_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrubs_SM_high$mean <- mean
newdat_shrubs_SM_high$lci <- lci
newdat_shrubs_SM_high$uci <- uci


#Plot
Shrub_seed_plot <- ggplot(newdat_shrubs_SM_high, aes(x = Seed_mass, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_shrubs_SM_low, aes(x = Seed_mass, y = mean), color = "blue") +
  geom_ribbon(data = newdat_shrubs_SM_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression(" ")) +
  xlab(expression("Seed mass (mg)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))


##### SLA #####
mean(Shrubs$SLA) #7.592492
sd(Shrubs$SLA) #2.498738
range(Shrubs$SLA) # 4.179176 12.149123

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR (870)
newdat_shrubs_SLA_low <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   SLA = seq(4, 13, by = 0.5), #range() 
                                   LDMC.std = 0, Height.std = 0, Seed_mass.std = 0, SR.std = -3.602942) %>% #predicting when other variables are at 0 which is their mean
  dplyr::mutate(SLA.std = (SLA - 7.592492)/2.498738) #back scaling
newdat_shrubs_SLA_low_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrubs_SLA_low, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_shrubs_SLA_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrubs_SLA_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrubs_SLA_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrubs_SLA_low$mean <- mean
newdat_shrubs_SLA_low$lci <- lci
newdat_shrubs_SLA_low$uci <- uci



##HIGH SR (980)
newdat_shrubs_SLA_high <- data.frame(Site = "average site", #predicting to the average site
                                    Species = "average species", #predicting to the average species
                                    SLA = seq(4, 13, by = 0.5), #range() 
                                    LDMC.std = 0, Height.std = 0, Seed_mass.std = 0, SR.std = 1.191129) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(SLA.std = (SLA - 7.592492)/2.498738) #back scaling)
#Predictions
newdat_shrubs_SLA_high_predictions<- plogis(posterior_epred(mTraitsShrubs, newdata = newdat_shrubs_SLA_high, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_shrubs_SLA_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_shrubs_SLA_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_shrubs_SLA_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_shrubs_SLA_high$mean <- mean
newdat_shrubs_SLA_high$lci <- lci
newdat_shrubs_SLA_high$uci <- uci


#Plot
Shrub_SLA_plot <- ggplot(newdat_shrubs_SLA_high, aes(x = SLA, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_shrubs_SLA_low, aes(x = SLA, y = mean), color = "blue") +
  geom_ribbon(data = newdat_shrubs_SLA_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression(" ")) +
  xlab(expression(SLA ~ (mm^2/mg^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))

Shrub_SLA_plot <- ggplot() +
  # High SLA
  geom_path(data = newdat_shrubs_SLA_high, aes(x = SLA, y = mean, color = "High")) +
  geom_ribbon(data = newdat_shrubs_SLA_high, aes(x = SLA, ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  # Low SLA
  geom_path(data = newdat_shrubs_SLA_low, aes(x = SLA, y = mean, color = "Low")) +
  geom_ribbon(data = newdat_shrubs_SLA_low, aes(x = SLA, ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression(" ")) +
  xlab(expression(SLA ~ (mm^2/mg^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) +
  scale_color_manual(values = c("High" = "red", "Low" = "blue")) +
  theme(legend.position = "right") +
  labs(color = "Solar radiation")


## Step 8) Plot ####
#Add dots in species into plots
shrub_species <- Shrubs %>% 
  select(Species, SLA, Seed_mass, LDMC, Height, SR)
#Get average values for all traits for each species
shrub_species <- shrub_species %>%
  group_by(Species) %>%
  summarise(SLA = mean(SLA), Seed_mass = mean(Seed_mass), LDMC = mean(LDMC), Height = mean(Height)) 
#Select for only Pimelea.alpine, Postanthera.cuneata, Acrothamnus.montanus, 
# Hovea.montana, and Podolobium.alpestre
shrub_species <- shrub_species %>%
  filter(Species %in% c("Olearia.phlogopappa.subsp..flavescens", 
                        "Grevillea.australis", 
                        "Euryomyrtus.ramosissima.subsp..ramosissima", 
                        "Acrothamnus.montanus", 
                        "Hovea.montana", 
                        "Podolobium.alpestre"))
unique(shrub_species$Species)

#Add in species
#SLA
Shrub_SLA_plot <- Shrub_SLA_plot +
  geom_point(data = shrub_species, aes(x = SLA, y = 0.25), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = shrub_species, aes(x = SLA, y = 0.25, label = Species), 
                           direction = "y", force = 5, size = 3,max.overlaps = 100)
#LDMC
Shrub_LDMC_plot <- Shrub_LDMC_plot +
  geom_point(data = shrub_species, aes(x = LDMC, y = 0.25), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = shrub_species, aes(x = LDMC, y = 0.25, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 100)
#Height
Shrub_height_plot <- Shrub_height_plot +
  geom_point(data = shrub_species, aes(x = Height, y = 0.25), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = shrub_species, aes(x = Height, y = 0.25, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 100)
#Seed mass
Shrub_seed_plot <- Shrub_seed_plot +
  geom_point(data = shrub_species, aes(x = Seed_mass, y = 0.25), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = shrub_species, aes(x = Seed_mass, y = 0.25, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 100)

#Add labels
TrShrub_fixed <- TrShrub_fixed+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
Shrub_rand <- Shrub_rand+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)

Shrub_LDMC_plot <- Shrub_LDMC_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
Shrub_SLA_plot <- Shrub_SLA_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)
Shrub_height_plot <- Shrub_height_plot+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)
Shrub_seed_plot <- Shrub_seed_plot+
  annotate("text", x = Inf, y = Inf, label = "(f)", vjust = 1, hjust = 1, size = 6)

#Plot
Shrub_trait_plot <- grid.arrange(TrShrub_fixed, 
                                 Shrub_rand,
                                 Shrub_LDMC_plot, 
                                 Shrub_SLA_plot,
                                 Shrub_height_plot, 
                                 Shrub_seed_plot,
                                 #layout_matrix = layout_matrix,
                                 ncol = 2)

ggsave("output/Traits/Bayesian/Shrubs/Shrub traits.png", plot = Shrub_trait_plot, width = 11, height = 12) 


########### GRAMINOID SPECIES ###################


glimpse(Graminoid)
hist(Graminoid$SR)
#n=19

## Step 1) Check correlation #####
#Only put variables in the model that have biological meaning and this goes the same for interactions
Graminoid %>%
  dplyr::select(-c(Code, Site, Quadrat, Origin, Species, Growth_form, qCover, Cover)) %>%
  GGally::ggpairs()
#low correlation <0.6 between traits


## Step 2) Standardise fixed and random variables #####
Graminoid$SR.std = as.vector(scale(Graminoid$SR))
Graminoid$SLA.std = as.vector(scale(Graminoid$SLA))
Graminoid$LDMC.std = as.vector(scale(Graminoid$LDMC))
Graminoid$Seed_mass.std = as.vector(scale(Graminoid$Seed_mass))
Graminoid$Height.std = as.vector(scale(Graminoid$Height))
glimpse(Graminoid)


## Step 3) Model #####
mTraitsGraminoid <- brms::brm(qCover ~
                                LDMC.std + SLA.std + Height.std + Seed_mass.std + SR.std +
                                LDMC.std:SR.std + Height.std:SR.std + Seed_mass.std:SR.std + SLA.std:SR.std +
                                (1|Site) + (1|Species),
                           family = gaussian,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           seed = 1234,
                           control = list(adapt_delta = 0.9, max_treedepth=15),
                           data = Graminoid)

## Step 4) Validate model ####
#Check model
summary(mTraitsGraminoid)
plot(mTraitsGraminoid)
pp_check(mTraitsGraminoid) 


## Step 5) Fixed effects coefficient plot ####
TrGram_fixed <- as.data.frame(fixef(mTraitsGraminoid))
TrGram_fixed2 <- rownames_to_column(TrGram_fixed, var = "Variable")
TrGram_fixed2$Variable <- factor(TrGram_fixed2$Variable,
                                 levels = c("LDMC.std:SR.std","SLA.std:SR.std","Height.std:SR.std","Seed_mass.std:SR.std", 
                                            "LDMC.std", "SLA.std", "Height.std","Seed_mass.std","SR.std",
                                            "Intercept"))

TrGram_fixed <- ggplot(TrGram_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")

## Step 6) Random effects #####

# 1. Plot estimates of the standard deviation and CI of random effects to see which one explains the most in the unexplained variation 
summary(mTraitsGraminoid)
Gram_ran_DF <- data.frame(
  variable = c("Site", "Species"),
  estimate = c( 0.41, 0.90),
  lci = c(0.26 ,0.65),
  uci = c(0.59, 1.40)
)

Gram_rand <- ggplot(Gram_ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

ggsave("output/Traits/Bayesian/Graminoids/Graminoid traits (random effect).png", 
       plot = Gram_rand, width = 10, height = 7.07) 


# 2. Site
TrGram_randSite <- as.data.frame(ranef(mTraitsGraminoid, groups = "Site"))
TrGram_randSite2 <- rownames_to_column(TrGram_randSite, var = "Variable")
TrShrub_randSite2 <- TrGram_randSite2[order(TrGram_randSite2$Site.Estimate.Intercept), ]

Gram_Site_rand <- ggplot(TrGram_randSite2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")

ggsave("output/Traits/Bayesian/Graminoids/Graminoid traits (Site random).png", 
       plot = Gram_Site_rand, width = 10, height = 7.07) 

# 3. Species
TrGram_randSpp <- as.data.frame(ranef(mTraitsGraminoid, groups = "Species"))
TrGram_randSpp2 <- rownames_to_column(TrGram_randSpp, var = "Variable")
TrGram_randSpp2 <- TrGram_randSpp2[order(TrGram_randSpp2$Species.Estimate.Intercept), ]

Gram_Species_rand <- ggplot(TrGram_randSpp2, aes(x = Species.Estimate.Intercept, y = reorder(Variable, Species.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Species.Q2.5.Intercept, xmax = Species.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Traits/Bayesian/Graminoids/Graminoids traits (Species random).png", 
       plot = Gram_Species_rand, width = 10, height = 7.07) 


## Step 7) Plot predictions #######
### LOW & HIGH SR
mean(Graminoid$SR) #461002.7
sd(Graminoid$SR) #21394.58
range(Graminoid$SR) 
##LOW SR (at 850W/m2 - back calcuated (850 -mean)/SD = -2.285559)
(393523.9 - 461002.7)/21394.58 = -3.154014
#HIGH SR (at 950W/m2 - back calcuated (850 -mean)/SD = -2.285559)
(486082.5 - 461002.7)/21394.58 = 1.17225

##### LDMC #####
mean(Graminoid$LDMC) #422.0755
sd(Graminoid$LDMC) #65.2567
range(Graminoid$LDMC) #326.2048 505.9174

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
##LOW SR
newdat_gram_ldmc_low <- data.frame(Site = "average site", #predicting to the average site
                              Species = "average species", #predicting to the average species
                              LDMC = seq(320, 510, by = 10), #range() of LDMC
                              Seed_mass.std = 0, SLA.std = 0, 
                              Height.std =0, SR.std = -3.154014) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 422.0755)/65.2567) #back scaling
#Predict
newdat_gram_ldmc_low_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_ldmc_low, 
                                                           allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_gram_ldmc_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_ldmc_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_ldmc_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_ldmc_low$mean <- mean
newdat_gram_ldmc_low$lci <- lci
newdat_gram_ldmc_low$uci <- uci


## HIGH SR
newdat_gram_ldmc_high <- data.frame(Site = "average site", #predicting to the average site
                               Species = "average species", #predicting to the average species
                               LDMC = seq(320, 510, by = 10), #range() of LDMC
                               Seed_mass.std = 0, SLA.std = 0, 
                               Height.std =0, SR.std = 1.17225) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 422.0755)/65.2567) #back scaling
#Predict
newdat_gram_ldmc_high_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_ldmc_high, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_gram_ldmc_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_ldmc_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_ldmc_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_ldmc_high$mean <- mean
newdat_gram_ldmc_high$lci <- lci
newdat_gram_ldmc_high$uci <- uci



### LDMC PLOT
Gram_LDMC_plot <- ggplot(newdat_gram_ldmc_high, aes(x = LDMC, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_gram_ldmc_low, aes(x = LDMC, y = mean), color = "blue") +
  geom_ribbon(data = newdat_gram_ldmc_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0,0.9) +
  ylab(expression("Graminoid species proportional cover")) +
  xlab(expression(LDMC ~ (mg/g^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))




##### Seed mass #####
mean(Graminoid$Seed_mass) #0.8335367
sd(Graminoid$Seed_mass) #1.597606
range(Graminoid$Seed_mass) #0.1266383 12.9700000

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat__gram_Seedmass_low <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   Seed_mass = seq(0.1, 13, by = 0.25), #range() of LDMC
                                   LDMC.std = 0, SLA.std = 0, 
                                   Height.std =0, SR.std = -3.154014) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = (Seed_mass - 0.8335367)/1.597606) #back scaling
#Model Predictions
newdat_gram_SM_low_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat__gram_Seedmass_low, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_gram_SM_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_SM_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_SM_low_predictions, 2, quantile, 0.975)
#Combine 
newdat__gram_Seedmass_low$mean <- mean
newdat__gram_Seedmass_low$lci <- lci
newdat__gram_Seedmass_low$uci <- uci


#HIGH SR
newdat_gram_Seedmass_high <- data.frame(Site = "average site", #predicting to the average site
                                    Species = "average species", #predicting to the average species
                                    Seed_mass = seq(0.1, 13, by = 0.25), #range() of LDMC
                                    LDMC.std = 0, SLA.std = 0, 
                                    Height.std =0, SR.std = 1.17225) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = (Seed_mass - 0.8335367)/1.597606) #back scaling
#Model predictions
newdat_gram_SM_high_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_Seedmass_high, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_gram_SM_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_SM_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_SM_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_Seedmass_high$mean <- mean
newdat_gram_Seedmass_high$lci <- lci
newdat_gram_Seedmass_high$uci <- uci


### Seed mass PLOT
Gram_SM_plot <- ggplot(newdat_gram_Seedmass_high, aes(x = Seed_mass, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat__gram_Seedmass_low, aes(x = Seed_mass, y = mean), color = "blue") +
  geom_ribbon(data = newdat__gram_Seedmass_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.9) +
  ylab(expression(" ")) +
  xlab(expression("Seed mass (mg)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))




##### SLA #####
mean(Graminoid$SLA) #12.2894
sd(Graminoid$SLA) #3.545692
range(Graminoid$SLA) # 4.375524 15.234369

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR (870)
newdat_gram_SLA_low <- data.frame(Site = "average site", #predicting to the average site
                                    Species = "average species", #predicting to the average species
                                    SLA = seq(4, 16, by = 0.5), #range() 
                                    LDMC.std = 0, Height.std = 0, 
                                  Seed_mass.std = 0, SR.std = -3.154014) %>% #predicting when other variables are at 0 which is their mean
  dplyr::mutate(SLA.std = (SLA - 12.2894)/3.545692) #back scaling
newdat_gram_SLA_low_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_SLA_low, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_gram_SLA_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_SLA_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_SLA_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_SLA_low$mean <- mean
newdat_gram_SLA_low$lci <- lci
newdat_gram_SLA_low$uci <- uci



##HIGH SR (980)
newdat_gram_SLA_high <- data.frame(Site = "average site", #predicting to the average site
                                     Species = "average species", #predicting to the average species
                                     SLA = seq(4, 16, by = 0.5), #range() 
                                     LDMC.std = 0, Height.std = 0, 
                                   Seed_mass.std = 0, SR.std = 1.17225) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(SLA.std = (SLA - 12.2894)/3.545692) #back scaling
#Predictions
newdat_gram_SLA_high_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_SLA_high, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_gram_SLA_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_SLA_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_SLA_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_SLA_high$mean <- mean
newdat_gram_SLA_high$lci <- lci
newdat_gram_SLA_high$uci <- uci


#Plot
Gram_SLA_plot <- ggplot(newdat_gram_SLA_high, aes(x = SLA, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_gram_SLA_low, aes(x = SLA, y = mean), color = "blue") +
  geom_ribbon(data = newdat_gram_SLA_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylab(expression(" ")) +
  ylim(0, 0.9) +
  xlab(expression(SLA ~ (mm^2/mg^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))




##### Height #####
mean(Graminoid$Height) #0.3929412
sd(Graminoid$Height) #0.2704131
range(Graminoid$Height) # 0.05 1.80

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR (870)
newdat_gram_height_low <- data.frame(Site = "average site", #predicting to the average site
                                  Species = "average species", #predicting to the average species
                                  Height = seq(0.05, 2, by = 0.05), #range() 
                                  LDMC.std = 0, SLA.std = 0, 
                                  Seed_mass.std = 0, SR.std = -3.154014) %>% #predicting when other variables are at 0 which is their mean
  dplyr::mutate(Height.std = (Height - 0.3929412)/0.2704131) #back scaling
newdat_gram_Height_low_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_height_low, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_gram_Height_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_Height_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_Height_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_height_low$mean <- mean
newdat_gram_height_low$lci <- lci
newdat_gram_height_low$uci <- uci



##HIGH SR (980)
newdat_gram_height_high <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   Height = seq(0.05, 2, by = 0.05), #range() 
                                   LDMC.std = 0, SLA.std = 0, 
                                   Seed_mass.std = 0, SR.std = 1.17225) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Height.std = (Height - 0.3929412)/0.2704131) #back scaling
#Predictions
newdat_gram_height_high_predictions<- plogis(posterior_epred(mTraitsGraminoid, newdata = newdat_gram_height_high, allow_new_levels = TRUE))
##Mean & CI
mean <- apply(newdat_gram_height_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_gram_height_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_gram_height_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_gram_height_high$mean <- mean
newdat_gram_height_high$lci <- lci
newdat_gram_height_high$uci <- uci


#Plot
Gram_height_plot <- ggplot(newdat_gram_height_high, aes(x = Height, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_gram_height_low, aes(x = Height, y = mean), color = "blue") +
  geom_ribbon(data = newdat_gram_height_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylab(expression("Graminoid species proportional cover")) +
  ylim(0, 0.9) +
  xlab(expression("Height (m)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))





## Step 8) Plot ####
### Step 8.2) Plot - add in species to plots
#Add dots in species into plots
gram_species <- Graminoid %>% 
  select(Species, SLA, Seed_mass, LDMC, Height, SR)
#Get average values for all traits for each species
gram_species <- gram_species %>%
  group_by(Species) %>%
  summarise(SLA = mean(SLA), Seed_mass = mean(Seed_mass), LDMC = mean(LDMC), Height = mean(Height)) 
#Select for only Pimelea.alpine, Postanthera.cuneata, Acrothamnus.montanus, 
# Hovea.montana, and Podolobium.alpestre
gram_species <- gram_species %>%
  filter(Species %in% c("Luzula.novae.cambriae", 
                         "Poa.fawcettiae", 
                         "Carex.breviculmis", 
                         "Anthosachne.scabra", 
                         "Rytidosperma.nudiflorum", 
                         "Rytidosperma.pallidum"))
unique(gram_species$Species)

#Add in species
#SLA
Gram_SLA_plot <- Gram_SLA_plot +
  geom_point(data = gram_species, aes(x = SLA, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = gram_species, aes(x = SLA, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3,max.overlaps = 10)
#LDMC
Gram_LDMC_plot <- Gram_LDMC_plot +
  geom_point(data = gram_species, aes(x = LDMC, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = gram_species, aes(x = LDMC, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)
#Height
Gram_height_plot <- Gram_height_plot +
  geom_point(data = gram_species, aes(x = Height, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = gram_species, aes(x = Height, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)
#Seed mass
Gram_SM_plot <- Gram_SM_plot +
  geom_point(data = gram_species, aes(x = Seed_mass, y = 0.25), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = gram_species, aes(x = Seed_mass, y = 0.25, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)

#Add labels
TrGram_fixed <- TrGram_fixed+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
Gram_rand <- Gram_rand+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)

Gram_LDMC_plot <- Gram_LDMC_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
Gram_SLA_plot <- Gram_SLA_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)
Gram_height_plot <- Gram_height_plot+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)
Gram_SM_plot <- Gram_SM_plot+
  annotate("text", x = Inf, y = Inf, label = "(f)", vjust = 1, hjust = 1, size = 6)

#Plot
Gram_trait_plot <- grid.arrange(TrGram_fixed, 
                                Gram_rand,
                                Gram_LDMC_plot, 
                                Gram_SLA_plot,
                                Gram_height_plot, 
                                Gram_SM_plot, 
                                ncol = 2)

ggsave("output/Traits/Bayesian/Graminoids/Graminoid traits.png", plot = Gram_trait_plot, width = 11, height = 12) 




############## FORB SPECIES #################


glimpse(Forbs)
hist(Forbs$SR)

## Step 1) Look for correlation between variables ############
#Only put variables in the model that have biological meaning and this goes the same for interactions
Forbs %>%
  dplyr::select(-c(Code, Site, Quadrat, Species, Origin, Growth_form, qCover)) %>%
  GGally::ggpairs()
#Low correlaton


## Step 2: Standardise fixed and random variables #########
#Scaled variables = (Shrubs$LDMC - mean)/sd(Shrubs$LDMC)
# I.e. Standardise variables (mean = 0 and SD = 1)

#Column for just scaled variables
Forbs$SR.std = as.vector(scale(Forbs$SR))
Forbs$SLA.std = as.vector(scale(Forbs$SLA))
Forbs$LDMC.std = as.vector(scale(Forbs$LDMC))
Forbs$Seed_mass.std = as.vector(scale(Forbs$Seed_mass))
Forbs$Height.std = as.vector(scale(Forbs$Height))
glimpse(Forbs)



## Step 2) Model ######## 

mTraitsForbs <- brms::brm(qCover ~
                            LDMC.std + SLA.std + Height.std + Seed_mass.std + SR.std +
                            LDMC.std:SR.std + Height.std:SR.std + Seed_mass.std:SR.std + SLA.std:SR.std +
                            (1|Site) + (1|Species),
                  family = gaussian,
                  iter = 4000,
                  cores = 4,
                  chains = 4,
                  seed = 1234,
                  control = list(adapt_delta = 0.9, max_treedepth=15),
                  data = Forbs)

## Step 3) Validate model  #########
#Check model
summary(mTraitsForbs)
plot(mTraitsForbs)
pp_check(mTraitsForbs) 


## Step 4) Fixed effects coefficient plot ########
mForb_fixed <- as.data.frame(fixef(mTraitsForbs))
mForb_fixed2 <- rownames_to_column(mForb_fixed, var = "Variable")
mForb_fixed2$Variable <- factor(mForb_fixed2$Variable,
                                levels = c("LDMC.std:SR.std","SLA.std:SR.std","Height.std:SR.std","Seed_mass.std:SR.std", 
                                           "LDMC.std", "SLA.std", "Height.std","Seed_mass.std","SR.std",
                                           "Intercept"))

custom_labels <- c("LDMC.std" = "LDMC", "SLA.std" = "SLA", "Height.std" = "Height", "Seed_mass.std" = "Seed mass", "SR.std" = "SR", 
                   "LDMC.std:SR.std" = "LDMC:SR", "SLA.std:SR.std" = "SLA:SR", "Height.std:SR.std" = "Height:SR", "Seed_mass.std:SR.std" = "Seed mass:SR",
                   "(Intercept)" = "Intercept")

Forb_fixed <- ggplot(mForb_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")


## Step 5) Random effects #######

# 1. Plot estimates of the standard deviation and CI of random effects to see which one explains the most in the unexplained variation 
summary(mTraitsForbs)
Forbs_ran_DF <- data.frame(
  variable = c("Site", "Species"),
  estimate = c(0.41, 0.80),
  lci = c(0.27 ,0.64),
  uci = c(0.61, 1.03)
)

Forbs_rand <- ggplot(Forbs_ran_DF, aes(x = estimate, y = variable)) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = lci, xmax = uci, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Estimated standard deviation")

ggsave("output/Traits/Bayesian/Forbs/Forbs traits (random effect).png", 
       plot = Gram_rand, width = 10, height = 7.07) 

# 2. Site effect
TrForb_randSite <- as.data.frame(ranef(mTraitsForbs, groups = "Site"))
TrForb_randSite2 <- rownames_to_column(TrForb_randSite, var = "Variable")
TrForb_randSite2 <- TrForb_randSite2[order(TrForb_randSite2$Site.Estimate.Intercept), ]

Forb_Site_rand <- ggplot(TrForb_randSite2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  #scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Logit coefficients")

ggsave("output/Traits/Bayesian/Forbs/Forbs traits (Site random).png", 
       plot = Forb_Site_rand, width = 10, height = 7.07) 


# 3. Species effect
TrForb_randSpp <- as.data.frame(ranef(mTraitsForbs, groups = "Species"))
TrForb_randSpp2 <- rownames_to_column(TrForb_randSpp, var = "Variable")
TrForb_randSpp2 <- TrForb_randSpp2[order(TrForb_randSpp2$Species.Estimate.Intercept), ]

Forb_Species_rand <- ggplot(TrForb_randSpp2, aes(x = Species.Estimate.Intercept, y = reorder(Variable, Species.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Species.Q2.5.Intercept, xmax = Species.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  #scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")


ggsave("output/Traits/Bayesian/Forbs/Forbs traits (Species random).png", 
       plot = Forb_Species_rand, width = 10, height = 7.07) 



## Step 6) Plot predictions ##########
### LOW & HIGH SR
mean(Forbs$SR) #462304.7
sd(Forbs$SR) #20020.31
range(Forbs$SR) 
##LOW SR (at 850W/m2 - back calcuated (850 -mean)/SD = -2.285559)
(393523.9 - 462304.7)/20020.31 = -3.435551
#HIGH SR (at 950W/m2 - back calcuated (850 -mean)/SD = -2.285559)
(486082.5 - 462304.7)/20020.31 = 1.187684



##### LDMC ######
mean(Forbs$LDMC) #210.9089
sd(Forbs$LDMC) #78.91859
range(Forbs$LDMC) #78.03423 613.63636

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_forbs_ldmc_low <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   LDMC = seq(78, 614, by = 5), #range() of LDMC
                                   Seed_mass.std = 0, SLA.std =0, 
                                   Height.std = 0, SR.std = -3.435551) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 210.9089)/78.91859)#back scaling
#Predictions
newdat_forb_ldmc_low_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forbs_ldmc_low, 
                                                         allow_new_levels = TRUE))
#Mean
mean <- apply(newdat_forb_ldmc_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_ldmc_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_ldmc_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_forbs_ldmc_low$mean <- mean
newdat_forbs_ldmc_low$lci <- lci
newdat_forbs_ldmc_low$uci <- uci


#HIGH SR
newdat_forbs_ldmc_HIGH <- data.frame(Site = "average site", #predicting to the average site
                                    Species = "average species", #predicting to the average species
                                    LDMC = seq(78, 614, by = 5), #range() of LDMC
                                    Seed_mass.std = 0, SLA.std = 0, Height.std = 0, SR.std = 1.187684) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(LDMC.std = (LDMC - 210.9089)/78.91859)#back scaling
#Predictions
newdat_forb_ldmc_high_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forbs_ldmc_HIGH, 
                                                          allow_new_levels = TRUE))
#Mean
mean <- apply(newdat_forb_ldmc_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_ldmc_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_ldmc_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_forbs_ldmc_HIGH$mean <- mean
newdat_forbs_ldmc_HIGH$lci <- lci
newdat_forbs_ldmc_HIGH$uci <- uci


### Plot
Forb_LDMC_plot <- ggplot(newdat_forbs_ldmc_HIGH, aes(x = LDMC, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_forbs_ldmc_low, aes(x = LDMC, y = mean), color = "blue") +
  geom_ribbon(data = newdat_forbs_ldmc_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.5) +
  ylab(expression("Forb species proportional cover")) +
  xlab(expression(LDMC ~ (mg/g^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))





##### SLA ######
mean(Forbs$SLA) #15.36089
sd(Forbs$SLA) #5.472301
range(Forbs$SLA) #2.223691 25.601682

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_forbs_sla_low <- data.frame(Site = "average site", #predicting to the average site
                             Species = "average species", #predicting to the average species
                             SLA = seq(2, 26, by = 0.25), #range() of LDMC
                             Seed_mass.std = 0, LDMC.std = 0, 
                             Height.std = 0, SR.std = -3.435551) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(SLA.std = ((SLA - 15.38108)/5.496649)) #back scaling
#Predictions
newdat_forb_SLA_low_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forbs_sla_low, 
                                                         allow_new_levels = TRUE))
#Mean
mean <- apply(newdat_forb_SLA_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_SLA_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_SLA_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_forbs_sla_low$mean <- mean
newdat_forbs_sla_low$lci <- lci
newdat_forbs_sla_low$uci <- uci


#HIGH SR
newdat_forbs_sla_HIGH <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   SLA = seq(2, 26, by = 0.25), #range() of LDMC
                                   Seed_mass.std = 0, LDMC.std = 0, 
                                   Height.std = 0, SR.std = 1.187684) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(SLA.std = ((SLA - 15.38108)/5.496649)) #back scaling
#Predictions
newdat_forb_SLA_high_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forbs_sla_HIGH, 
                                                         allow_new_levels = TRUE))
#Mean
mean <- apply(newdat_forb_SLA_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_SLA_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_SLA_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_forbs_sla_HIGH$mean <- mean
newdat_forbs_sla_HIGH$lci <- lci
newdat_forbs_sla_HIGH$uci <- uci



### Plot
Forb_SLA_plot <- ggplot(newdat_forbs_sla_HIGH, aes(x = SLA, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_forbs_sla_low, aes(x = SLA, y = mean), color = "blue") +
  geom_ribbon(data = newdat_forbs_sla_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.5) +
  ylab(expression(" ")) +
  xlab(expression(SLA ~ (mm^2/mg^-1))) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))





##### Height ####
mean(Forbs$Height) #0.4377948
sd(Forbs$Height) #0.2526176
range(Forbs$Height) #0.019 1.300

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_forb_H_height_low <- data.frame(Site = "average site", #predicting to the average site
                                  Species = "average species", #predicting to the average species
                                  Height = seq(0.015, 1.4, by = 0.5), #range() of LDMC
                                  Seed_mass.std = 0,LDMC.std = 0, 
                                  SLA.std = 0, SR.std = -3.435551) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Height.std = ((Height - 0.4377948)/0.2526176)) #back transforming and unscaling
#Predict
newdat_forb_height_low_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forb_H_height_low, 
                                                          allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_forb_height_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_height_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_height_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_forb_H_height_low$mean <- mean
newdat_forb_H_height_low$lci <- lci
newdat_forb_H_height_low$uci <- uci


#HIGH SR
newdat_forb_H_height_high <- data.frame(Site = "average site", #predicting to the average site
                                   Species = "average species", #predicting to the average species
                                   Height = seq(0.015, 1.4, by = 0.5), #range() of LDMC
                                   Seed_mass.std = 0, LDMC.std = 0, 
                                   SLA.std = 0, SR.std = 1.187684) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Height.std = ((Height - 0.4377948)/0.2526176)) 
#Predict
newdat_forb_height_high_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forb_H_height_high, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_forb_height_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_height_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_height_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_forb_H_height_high$mean <- mean
newdat_forb_H_height_high$lci <- lci
newdat_forb_H_height_high$uci <- uci


### Height PLOT
Forb_height_plot <- ggplot(newdat_forb_H_height_high, aes(x = Height, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_forb_H_height_low, aes(x = Height, y = mean), color = "blue") +
  geom_ribbon(data = newdat_forb_H_height_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.5) +
  ylab(expression("Forb species proportional cover")) +
  xlab(expression("Height (m)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))




#### Seed mass ####
mean(Forbs$Seed_mass) #1.768466
sd(Forbs$Seed_mass) #1.764379
range(Forbs$Seed_mass) #0.007431352 5.409800000

#Create new dataset to predict to the average site when other variables at at 0 (their mean)
#LOW SR
newdat_forb_Seed_mass_low <- data.frame(Site = "average site", #predicting to the average site
                                     Species = "average species", #predicting to the average species
                                     Seed_mass = seq(0.006, 6, by = 0.05), #range() of LDMC
                                     Height.std = 0, LDMC.std = 0, 
                                     SLA.std = 0, SR.std = -3.435551) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = ((Seed_mass - 1.766151)/1.764351)) #back transforming and unscaling
#Model predictions
newdat_forb_SM_low_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_forb_Seed_mass_low, 
                                                          allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_forb_SM_low_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_SM_low_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_SM_low_predictions, 2, quantile, 0.975)
#Combine 
newdat_forb_Seed_mass_low$mean <- mean
newdat_forb_Seed_mass_low$lci <- lci
newdat_forb_Seed_mass_low$uci <- uci



#HIGH SR
newdat_Forbs_Seed_mass_high <- data.frame(Site = "average site", #predicting to the average site
                                      Species = "average species", #predicting to the average species
                                      Seed_mass = seq(0.006, 6, by = 0.05), #range() of LDMC
                                      Height.std = 0, LDMC.std = 0, 
                                      SLA.std = 0, SR.std = 1.187684) %>% #predicting LDMC when SLA, height and SR are at 0 which is their mean
  dplyr::mutate(Seed_mass.std = ((Seed_mass - 1.766151)/1.764351)) #back transforming and unscaling
#Model predictions
newdat_forb_SM_high_predictions<- plogis(posterior_epred(mTraitsForbs, newdata = newdat_Forbs_Seed_mass_high, allow_new_levels = TRUE))
#Calcuate mean and CI
mean <- apply(newdat_forb_SM_high_predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(newdat_forb_SM_high_predictions, 2, quantile, 0.025)
uci <- apply(newdat_forb_SM_high_predictions, 2, quantile, 0.975)
#Combine 
newdat_Forbs_Seed_mass_high$mean <- mean
newdat_Forbs_Seed_mass_high$lci <- lci
newdat_Forbs_Seed_mass_high$uci <- uci



### Seed mass PLOT
Forb_SM_plot <- ggplot(newdat_Forbs_Seed_mass_high, aes(x = Seed_mass, y = mean)) +
  geom_path(color = "red") +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "red") +
  geom_path(data = newdat_forb_Seed_mass_low, aes(x = Seed_mass, y = mean), color = "blue") +
  geom_ribbon(data = newdat_forb_Seed_mass_low, aes(ymin = lci, ymax = uci), alpha = 0.2, fill = "blue") +
  theme_cowplot() +
  ylim(0, 0.5) +
  ylab(expression(" ")) +
  xlab(expression("Seed mass (mg)")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))


## Step 8) Plot ####
### Step 8.2) Plot - add in species to plots
#Add dots in species into plots
forb_species <- Forbs %>% 
  select(Species, SLA, Seed_mass, LDMC, Height, SR)
#Get average values for all traits for each species
forb_species <- forb_species %>%
  group_by(Species) %>%
  summarise(SLA = mean(SLA), Seed_mass = mean(Seed_mass), LDMC = mean(LDMC), Height = mean(Height)) 
#Select for only Pimelea.alpine, Postanthera.cuneata, Acrothamnus.montanus, 
# Hovea.montana, and Podolobium.alpestre
forb_species <- forb_species %>%
  filter(Species %in% c("Hypochaeris.radicata", 
                        "Celmisia.costiniana", 
                        "Craspedia.aurantia", 
                        "Ewartia.nubigena", 
                        "Stelleria.pungens", 
                        "Leptorhynchos.squamatus.subsp..alpinus"))
unique(forb_species$Species)

#Add in species
#SLA
Forb_SLA_plot <- Forb_SLA_plot +
  geom_point(data = forb_species, aes(x = SLA, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = forb_species, aes(x = SLA, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3,max.overlaps = 10)
#LDMC
Forb_LDMC_plot <- Forb_LDMC_plot +
  geom_point(data = forb_species, aes(x = LDMC, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = forb_species, aes(x = LDMC, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)
#Height
Forb_height_plot <- Forb_height_plot +
  geom_point(data = forb_species, aes(x = Height, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = forb_species, aes(x = Height, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)
#Seed mass
Forb_SM_plot <- Forb_SM_plot +
  geom_point(data = forb_species, aes(x = Seed_mass, y = 0.10), 
             size = 3, shape = 21, fill = "black") +
  ggrepel::geom_text_repel(data = forb_species, aes(x = Seed_mass, y = 0.10, label = Species), 
                           direction = "y", force = 5, size = 3, max.overlaps = 10)

#Add labels
Forb_fixed <- Forb_fixed+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
Forbs_rand <- Forbs_rand+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)

Forb_LDMC_plot <- Forb_LDMC_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
Forb_SLA_plot <- Forb_SLA_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)
Forb_height_plot <- Forb_height_plot+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)
Forb_SM_plot <- Forb_SM_plot+
  annotate("text", x = Inf, y = Inf, label = "(f)", vjust = 1, hjust = 1, size = 6)



# Plot
plot_forb <- grid.arrange(Forb_fixed, 
                          Forbs_rand,
                          Forb_LDMC_plot, 
                          Forb_SLA_plot, 
                          Forb_height_plot, 
                          Forb_SM_plot,
                          ncol = 2)

ggsave("output/Traits/Bayesian/Forbs/Forbs traits.png", 
       plot = plot_forb,  width = 11, height = 12)




## APPENDIX ########################################################################
##### Plot predictions vs observed for appendix 
#Shrubs
pdf("output/Appendix/Shrub traits model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mTraitsShrubs))
plot(preds$Estimate ~ mTraitsShrubs$data$qCover,
     xlab = "Shrub observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')

#Graminoids
pdf("output/Appendix/Graminoid traits model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mTraitsGraminoid))
plot(preds$Estimate ~ mTraitsGraminoid$data$qCover,
     xlab = "Graminoid observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')

#Forbs
pdf("output/Appendix/Forb traits model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mTraitsForbs))
plot(preds$Estimate ~ mTraitsForbs$data$qCover,
     xlab = "Forb observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')



