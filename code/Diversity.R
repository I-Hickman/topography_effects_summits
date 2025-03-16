##### DIVERSITY GLMM #######
library(performance)
library(see)
library(adespatial)
library(datawizard)
library(MASS)
library(report) 
library(ggeffects)
library(cowplot)
library(sjPlot)
library(ggeffects)
library(glmmTMB)
library(webshot)
library(vegan)
library(broom.mixed)
library(jtools)
library(effects)
library(gridExtra)
library(ggplot2)
library(brms)
library(dplyr)
library(tibble)


###
##
#
# CONTENTS
# 1. Load data and scale and center variables
# 2. Richness vs solar radiation
# 3. Diversity Hill 1 vs solar radiation
# 4. Diversity Hill 2 vs Solar radiation
# 5. Save all graphs
###
##
#


### 1) LOAD DATA AND SCALE VARIABLES #####

div_df <- read.csv("data/Diversity dataset.csv")
glimpse(div_df) 

div_df$SR.std = as.vector(scale(div_df$avg.SR)) #scales and centers variables


############ 1) RICHNESS ###########
#Look at the relationship between species richness and solar radiation
#Step 1) Plot the data
ggplot(div_df, aes(x = avg.SR, y = species_richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Species richness vs Solar radiation",
       x = "Solar radiation (standardised)",
       y = "Species richness") +
  theme_minimal()

#### Step 2) Build ecologiclly sound models that attemt to answer your question

mD0 <- brms::brm(species_richness ~ SR.std + (1 | Site),
                  family = poisson,
                  iter = 4000,
                  cores = 4,
                  chains = 4,
                  control = list(adapt_delta = 0.8, max_treedepth=10),
                  data = div_df)

#### Step 3) Validate model 
summary(mD0)
plot(mD0)
pp_check(mD0) 
prior_summary(mD0)


##### Step 4) Fixed effects coefficient plot
mRich_fixed <- as.data.frame(fixef(mD0))
mRich_fixed2 <- rownames_to_column(mRich_fixed, var = "Variable")
mRich_fixed2$Variable <- factor(mRich_fixed2$Variable,
                               levels = c("SR.std", "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "Intercept" = "Intercept")

Rich_fixed_plot <- ggplot(mRich_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")


#### Step 5) Random effects

mRich_rand <- as.data.frame(ranef(mD0))
mRich_rand2 <- rownames_to_column(mRich_rand, var = "Variable")
mRich_rand2 <- mRich_rand2[order(mRich_rand2$Site.Estimate.Intercept), ]

Rich_rand <- ggplot(mRich_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Diversity indices/Bayesian/D0/D0 (random).png", 
       plot = Rich_rand, width = 5, height = 3.5) 



#### Step 4) Plot predictions
mean(div_df$avg.SR) #462020.3
sd(div_df$avg.SR) #20615.8
range(div_df$avg.SR) #393523.9 486082.5

new_rich <- data.frame(Site = "average site", 
                       avg.SR = seq(393520, 486080, by = 1000)) %>% 
  dplyr::mutate(SR.std = (avg.SR - 462020.3)/20615.8)
#Predictions
predictions<- posterior_epred(mD0, newdata = new_rich,  
                              allow_new_levels = TRUE) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_rich$mean <- mean
new_rich$lci <- lci
new_rich$uci <- uci

### Plot
D0_plot <- ggplot(new_rich, aes(x = avg.SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = div_df, aes(x = avg.SR, y = species_richness), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylab(expression("D" ^ "0"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

#Add labels
Rich_fixed_plot <- Rich_fixed_plot+
  annotate("text", x = Inf, y = Inf, label = "(a)", vjust = 1, hjust = 1, size = 6)
D0_plot <- D0_plot+
  annotate("text", x = Inf, y = Inf, label = "(b)", vjust = 1, hjust = 1, size = 6)


#PLot all togeher
plot_D0 <- grid.arrange(Rich_fixed_plot, D0_plot, ncol = 2)

ggsave("output/Diversity indices/Bayesian/D0/D0 vs solar radiation.png", 
       plot = plot_D0, width = 10, height = 3.5)




######## 3) HILL 1 #########
#look at relationship between solar radiation and Hill1 using ggplot
ggplot(div_df, aes(x = avg.SR, y = Hill1)) + 
  geom_smooth(method = "lm") +
  geom_point() +
  theme_cowplot() +
  ylab(expression("Hill" ^ "1"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12))

#### Step 1) Build ecologically sound model
mD1 <- brms::brm(Hill1 ~ SR.std + (1 | Site), 
                      family = Gamma(link = "log"),
                      iter = 4000,
                      cores = 4,
                      chains = 4,
                      control = list(adapt_delta = 0.8, max_treedepth=10),
                      data = div_df)

#Check model
summary(mD1)
plot(mD1)
pp_check(mD1)
prior_summary(mD1)



##### Fixed effects coefficient plot
D1_fixed <- as.data.frame(fixef(mD1))
D1_fixed2 <- rownames_to_column(D1_fixed, var = "Variable")
D1_fixed2$Variable <- factor(D1_fixed2$Variable,
                              levels = c("SR.std", "Intercept"))

#Plot
custom_labels <- c("SR.std" = "SR",
                   "Intercept" = "Intercept")

D1_fixed_plot <- ggplot(D1_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")



#### Random effects coefficient plot
D1_rand <- as.data.frame(ranef(mD1))
D1_rand2 <- rownames_to_column(D1_rand, var = "Variable")
D1_rand2 <- D1_rand2[order(D1_rand2$Site.Estimate.Intercept), ]

D1_rand <- ggplot(D1_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")

ggsave("output/Diversity indices/Bayesian/D1/D1 (random).png", 
       plot = D1_rand, width = 5, height = 3.5) 


#### Step 3) Plot predictions from model
#new dataframe
new_hill1 <- data.frame(Site = "average site",
                        avg.SR = seq(393520, 486080, by = 1000)) %>% 
  dplyr::mutate(SR.std = (avg.SR - 462020.3)/20615.8)
#Predictions
predictions<- posterior_epred(mD1, newdata = new_hill1,  
                              allow_new_levels = TRUE) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_hill1$mean <- mean
new_hill1$lci <- lci
new_hill1$uci <- uci

### Plot
D1_plot <- ggplot(new_hill1, aes(x = avg.SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = div_df, aes(x = avg.SR, y = Hill1), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylab(expression("D" ^ "1"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

#Add labels
D1_fixed_plot <- D1_fixed_plot+
  annotate("text", x = Inf, y = Inf, label = "(c)", vjust = 1, hjust = 1, size = 6)
D1_plot <- D1_plot+
  annotate("text", x = Inf, y = Inf, label = "(d)", vjust = 1, hjust = 1, size = 6)

#PLot all togeher
plot_D1 <- grid.arrange(D1_fixed_plot, D1_plot, ncol = 2)

ggsave("output/Diversity indices/Bayesian/D1/D1 vs solar radiation.png", 
       plot = plot_D1, width = 10, height = 3.5)




########## 4) HILL 2 ###########

#### Step 1) Build ecologically sound model

mD2 <- brms::brm(Hill2 ~ SR.std + (1 | Site),
                  family = Gamma(link = "log"),
                  iter = 4000,
                  cores = 4,
                  chains = 4,
                  control = list(adapt_delta = 0.8, max_treedepth=10),
                  data = div_df)

#### Step 2) Validate model 
summary(mD2)
plot(mD2)
pp_check(mD2) 


##### Step 4) Fixed effects coefficient plot
mD2_fixed <- as.data.frame(fixef(mD2))
mD2_fixed2 <- rownames_to_column(mD2_fixed, var = "Variable")
mD2_fixed2$Variable <- factor(mD2_fixed2$Variable,
                               levels = c("SR.std", "Intercept"))

custom_labels <- c("SR.std" = "SR",
                   "Intercept" = "Intercept")

D2_fixed <- ggplot(mD2_fixed2, aes(x = Estimate, y = Variable)) +
  geom_point(shape = 16, size = 2.5) +         # Point for the estimate
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5, height = 0)) +  # Confidence intervals
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    # Add a dashed line at x = 0
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")


#### Step 5) Random effects

mD2_rand <- as.data.frame(ranef(mD2))
mD2_rand2 <- rownames_to_column(mD2_rand, var = "Variable")
mD2_rand2 <- mD2_rand2[order(mD2_rand2$Site.Estimate.Intercept), ]

D2_rand <- ggplot(mD2_rand2, aes(x = Site.Estimate.Intercept, y = reorder(Variable, Site.Estimate.Intercept))) +
  geom_point(shape = 16, size = 2.5) +        
  geom_errorbarh(aes(xmin = Site.Q2.5.Intercept, xmax = Site.Q97.5.Intercept, height = 0)) +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +    
  scale_y_discrete(labels = custom_labels) +
  scale_x_continuous() +
  theme_half_open() +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 12)) +
  xlab("Coefficients")


ggsave("output/Diversity indices/Bayesian/D2/D2 (random).png", 
       plot = D2_rand, width = 5, height = 3.5) 


#### Step 4) Plot predictions

mean(div_df$avg.SR) #462020.3
sd(div_df$avg.SR) #20615.8
range(div_df$avg.SR) #393523.9 486082.5
#new dataframe
new_D2 <- data.frame(Site = "average site", 
                     avg.SR = seq(393520, 486080, by = 1000)) %>% 
  dplyr::mutate(SR.std = (avg.SR - 462020.3)/20615.8)
#Predictions
predictions<- posterior_epred(mD2, newdata = new_D2,  
                              allow_new_levels = TRUE) 
#Mean
mean <- apply(predictions, 2, mean)
#Confidence limits from predictions (the 95th credible intervals)
lci <- apply(predictions, 2, quantile, 0.025)
uci <- apply(predictions, 2, quantile, 0.975)
#Combine 
new_D2$mean <- mean
new_D2$lci <- lci
new_D2$uci <- uci

### Plot
D2_plot <- ggplot(new_D2, aes(x = avg.SR, y = mean)) + 
  geom_path() +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_point(data = div_df, aes(x = avg.SR, y = Hill2), 
             shape = 21, fill = "lightgrey", alpha = 0.2, colour = "black")+
  geom_line(size = 0.75) +
  theme_cowplot()+
  ylab(expression("D" ^ "2"))+
  xlab(expression("Solar radiation W/m" ^ "2")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12)) 

#Add labels
D2_fixed <- D2_fixed+
  annotate("text", x = Inf, y = Inf, label = "(e)", vjust = 1, hjust = 1, size = 6)
D2_plot <- D2_plot+
  annotate("text", x = Inf, y = Inf, label = "(f)", vjust = 1, hjust = 1, size = 6)

#PLot all togeher
plot_D2 <- grid.arrange(D2_fixed, D2_plot, ncol = 2)

ggsave("output/Diversity indices/Bayesian/D2/D2 vs solar radiation.png", 
       plot = plot_D2, width = 10, height = 3.5)


#####SAVE ALL ########
plot_all <- grid.arrange(plot_D0, plot_D1, plot_D2, nrow = 3, ncol = 1)

ggsave("output/Diversity indices/Bayesian/Diversity indices vs solar radiation.png", 
       plot = plot_all, width = 11, height = 12)


##### APPENDIX: Plot predictions vs observed for appendix #####
#D0
pdf("output/Appendix/mD0 model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mD0))
plot(preds$Estimate ~ mD0$data$species_richness,  
     xlab = "D0 observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 

#D1
pdf("output/Appendix/mD1 model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mD1))
plot(preds$Estimate ~ mD1$data$Hill1,  
     xlab = "D1 observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 
#D2
pdf("output/Appendix/mD2 model predictions.pdf", width = 5, height = 3.5)
preds <- as.data.frame(fitted(mD2))
plot(preds$Estimate ~ mD2$data$Hill2,  
     xlab = "D2 observations",
     ylab = "Model predictions")
abline(0, 1, col= 'red')
dev.off() 




