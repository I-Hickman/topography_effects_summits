#Calculating Diversity with hill numnbers
#Vegan package and use divesrity function then transform to hill numbers

library(vegan)
library(dplyr)
library(tidyr)

##### Step 1: Calculating alpha diversity indices
# Step 1: put together all the variables to form a code and then calcuaute hills and then put onto richness DF

DF_2022 <- read.csv("data/2022 Summit floristic.csv")

#Step 2: Creating matrix 

DF_test <- DF_2022[,c(1,14:122)] #gather your code and species columns
trans_test <- t(DF_test) #transpose it so that species runs down and code runs along the top
colnames(trans_test) <-  trans_test[1,] #make row 1 your top row
trans_test2 <-  trans_test[-1,] #remove the top row that is wrong
trans_test3 <- apply(trans_test2, 2, as.numeric) #change dataset to numerical

# Step 3: Calculating Hill's numbers

Shannon_div <- vegan::diversity(trans_test3, index = "shannon", equalize.groups = FALSE,   
                                MARGIN = 2, base = exp(1)) # Margin 2 assuming you have done the t() function on your raw % data
hill_1 <- exp(Shannon_div) # this is your Hill's number order = 1 vector for every site. 

simps <- vegan::diversity(trans_test3, index = "invsimpson", equalize.groups = FALSE,   
                          MARGIN = 2, base = exp(1)) 
hill_2 <- 1/(simps) # this is your Hill # 2 vector 

#Turn output into dataframes
Hill_1 <- as.data.frame(hill_1)
Hill_2 <- as.data.frame(hill_2)

# Change the rows into column 1
Hill_1 <- cbind(Code = rownames(Hill_1), Hill_1)
rownames(Hill_1) <- NULL
Hill_2 <- cbind(Code = rownames(Hill_2), Hill_2)
rownames(Hill_2) <- NULL


#Join onto the species richness data
Full_diversity_DF3 <-read.csv("data/Richness dataset.csv")
Full_diversity_DF4 <- Full_diversity_DF3 %>% 
  left_join(Hill_1, by = "Code")
Full_diversity_DF4 <- Full_diversity_DF4 %>% 
  left_join(Hill_2, by = "Code")

#Save dataset for future use
write.csv(Full_diversity_DF4, "data/Diversity dataset.csv", row.names = FALSE)
