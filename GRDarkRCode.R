############# Growth Rate Dark data only ####################


##I will load my necessary libraries...

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)

#loading in the data, that might be useful
GRData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/MScDiatomRCode/GRData.csv", header = TRUE)

names(GRData)[1] <- "ID"


#### Now I want to subset the data to only the dark data:

# Now I want to subset the full dataframe to only include the dark treatment

GRData_dark <- GRData %>% filter(treatment == "Dark")

head(GRData_dark)

## Ok, I think that worled nicely



GRBDarkMod1 <- lmer(growthrateB ~ temperature*day + (1|strain) + (1|ID), data = GRData_dark)
summary(GRBDarkMod1)


plot(fitted(GRBDarkMod1), residuals(GRBDarkMod1))
abline(h = 0, col = "red")

qqnorm(residuals(GRBDarkMod1))
qqline(residuals(GRBDarkMod1), col = "red")

# Deviance and AIC
deviance(GRBDarkMod1)
#-323.966
AIC(GRBDarkMod1)
#-309.966

# Overdispersion check
overdispersion <- deviance(GRBDarkMod1) / df.residual(GRBDarkMod1)
overdispersion
#-2.866956


# Cook's distance
cooksd <- cooks.distance(GRBDarkMod1)
plot(cooksd, type = "h")
abline(h = 4 / nrow(GRData_dark), col = "red") # Rule of thumb for identifying influential points
#??




###################################################################################


GRData_light <- GRData %>% filter(treatment == "Control")


head(GRData_light)

## Ok, I think that worled nicely



GRBLightMod1 <- lmer(growthrateB ~ temperature*day + (1|strain) + (1|ID), data = GRData_light)
summary(GRBLightMod1)


plot(fitted(GRBLightMod1), residuals(GRBLightMod1))
abline(h = 0, col = "red")

qqnorm(residuals(GRBLightMod1))
qqline(residuals(GRBLightMod1), col = "red")

# Deviance and AIC
deviance(GRBLightMod1)
#-219.3406
AIC(GRBLightMod1)
#-205.3406

# Overdispersion check
overdispersion <- deviance(GRBLightMod1) / df.residual(GRBLightMod1)
overdispersion
#-2.642658






























