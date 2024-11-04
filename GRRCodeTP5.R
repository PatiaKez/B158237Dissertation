############# FUll Growth Rate R Code with TP5 (Finally) ##################

##I will load my necessary libraries...

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)

#loading in the data, that might be useful
GRData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/B158237Dissertation/GRData5.csv", header = TRUE)

names(GRData)[1] <- "ID"


##Very first plot
plot(growthrateA~day, data =GRData)

plot(growthrateB~day, data = GRData)

### Now I will properly visualise the data:

## Raw Growth Rates

ggplot(GRData, aes(x = day, y = growthrateA, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Growth rate after re-inoculation into Light",
       x = "Timepoint (days)",
       y = "Growth Rate") +
  theme_minimal()

## Corrected Growth Rates

ggplot(GRData, aes(x = day, y = growthrateB, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Growth rate after re-inoculation into Light",
       x = "Timepoint (days)",
       y = "Growth Rate") +
  theme_minimal()


### and on to stats:

library(lmerTest)

#ok, let's give this a shot

GRAMod1 <- lm(growthrateA ~ strain + temperature + treatment + day, data = GRData)
summary(GRAMod1)

#let's take a look at model diagnostics
plot(GRAMod1)

#hmmm, doesn't seem too bad, but let's check the data distribution:

hist(GRData$growthrateA)

#Oh SHIT! That's normally distributed!

#ok, twist, wasn't expecting that!

#ok, well then...let's look at growthrateB

hist(GRData$growthrateB)

###############

#well, let's try a more sophisticated model then:

GRAMod2 <- lm(growthrateA ~ temperature*treatment*day + strain, data = GRData)
summary(GRAMod2)


GRAMod3 <- lmer(growthrateA ~ temperature*treatment*day + (1|strain) + (1|ID), data = GRData)
summary(GRAMod3)
## wow, this model works...
# was not expecting it

#well, slow down, let's look at plot diagnostics:

plot(GRAMod3)

plot(fitted(GRAMod3), residuals(GRAMod3))
abline(h = 0, col = "red")

qqnorm(residuals(GRAMod3))
qqline(residuals(GRAMod3), col = "red")

# Deviance and AIC
deviance(GRAMod3)
AIC(GRAMod3)

# Overdispersion check
overdispersion <- deviance(GRAMod3) / df.residual(GRAMod3)
overdispersion
#hmm... could be too overdispersed...


# Cook's distance
cooksd <- cooks.distance(GRAMod3)
plot(cooksd, type = "h")
abline(h = 4 / nrow(GRData), col = "red") # Rule of thumb for identifying influential points
#??


# Install and load the car package if not already installed
install.packages("car")
library(car)

# Calculate VIF
vif(GRData, type = "predictor")


#######################

### Let's repeat for the Corrected growth rates:

GRBMod1 <- lm(growthrateB ~ strain + temperature + treatment + day, data = GRData)
summary(GRBMod1)

#Checking the data distribution:
hist(GRData$growthrateB)
#NORMAL!!!

plot(fitted(GRBMod1), residuals(GRBMod1))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod1))
qqline(residuals(GRBMod1), col = "red")

# Deviance and AIC
deviance(GRBMod1)
#-564.375
AIC(GRBMod1)
#-542.375



GRBMod3 <- lmer(growthrateB ~ temperature*treatment*day + (1|strain) + (1|ID), data = GRData)
summary(GRBMod3)


plot(fitted(GRBMod3), residuals(GRBMod3))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod3))
qqline(residuals(GRBMod3), col = "red")

# Deviance and AIC
deviance(GRBMod3)
#-564.375
AIC(GRBMod3)
#-542.375

# Overdispersion check
overdispersion <- deviance(GRBMod3) / df.residual(GRBMod3)
overdispersion
#-2.46452


# Cook's distance
cooksd <- cooks.distance(GRBMod3)
plot(cooksd, type = "h")
abline(h = 4 / nrow(GRData), col = "red") # Rule of thumb for identifying influential points
#??

### This is so beautiful I could cry, I am crying 

## Just trying this simplyfied version?
GRBMod3b <- lmer(growthrateB ~ treatment*day + temperature + (1|strain) + (1|ID), data = GRData)
summary(GRBMod3b)


plot(fitted(GRBMod3b), residuals(GRBMod3b))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod3b))
qqline(residuals(GRBMod3b), col = "red")

# Deviance and AIC
deviance(GRBMod3b)
#-585.5207
AIC(GRBMod3b)
#-569.5207

# Overdispersion check
overdispersion <- deviance(GRBMod3b) / df.residual(GRBMod3b)
overdispersion
#-2.523796

## GRBMod3b does not fit quite as well as GRBMod3


GRBMod3c <- lmer(growthrateB ~ temperature*day + treatment + (1|strain) + (1|ID), data = GRData)
summary(GRBMod3c)

plot(fitted(GRBMod3c), residuals(GRBMod3c))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod3c))
qqline(residuals(GRBMod3c), col = "red")

# Deviance and AIC
deviance(GRBMod3c)
#-478.1676
AIC(GRBMod3c)
#-462.1676

# Overdispersion check
overdispersion <- deviance(GRBMod3c) / df.residual(GRBMod3c)
overdispersion
#-2.061067

## ok fuck... this model seems to fit best...


####################################################################################


### Let's do a manner of post-hov test, including strain as a fixed effect:

GRBMod4 <- lmer(growthrateB ~ temperature*treatment*day + strain + (1|ID), data = GRData)
summary(GRBMod4)


plot(fitted(GRBMod4), residuals(GRBMod4))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod4))
qqline(residuals(GRBMod4), col = "red")

# Deviance and AIC
deviance(GRBMod4)
#-549.0447
AIC(GRBMod4)
#-521.0447

# Overdispersion check
overdispersion <- deviance(GRBMod4) / df.residual(GRBMod4)
overdispersion
#hmm... could be too overdispersed...


# Cook's distance
cooksd <- cooks.distance(GRBMod4)
plot(cooksd, type = "h")
abline(h = 4 / nrow(GR3Data), col = "red") 
# Rule of thumb for identifying influential points
#??






