### EB R Code -  Finding the best model ####

### ok, well biology, statistical rules, and logic dictate that I need to inlcude strain
#as a random effect

# furthermore it seems I might have to include replicate ID as a random effect, 

#So here I am assessing which model is the best fit for my data:

## First I need to inlcude the data loading, library loading, and some prep code:

##I will load my necessary libraries...

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)

#loading in the data, that might be useful
EBData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/B158237Dissertation/EBData.csv", header = TRUE)

names(EBData)[1] <- "ID"

## First I will create a new column and transform my percentage into proportion dead:

EBData$prop_dead <- EBData$perc_dead/100

## I want to convert temperature into an ordinal variable:

#EBData$temperature <- factor(EBData$temperature, levels = c(1, 4), ordered = TRUE)


### I also want to visulaise the data:

# Create the plot
ggplot(EBData, aes(x = day, y = prop_dead, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Proportion of Dead Cells Over Time",
       x = "Timepoint (weeks)",
       y = "Proportion of Dead Cells") +
  theme_minimal()

## lovely


#### ok, now i will start with the model I know works:

EBMod5 <- glm(cbind(dead,live) ~ temperature*treatment*day + strain, family = binomial(link="logit"), data = EBData)
summary(EBMod5)

plot(EBMod5)


plot(fitted(EBMod5), residuals(EBMod5))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod5))
qqline(residuals(EBMod5), col = "red")

# Deviance and AIC
deviance(EBMod5)
AIC(EBMod5)

# Overdispersion check
overdispersion <- deviance(EBMod5) / df.residual(EBMod5)
overdispersion
#hmm... could be too overdispersed...


#######################

##### Now for the model with strain as a random effect:
#Note, I named this Mod10 so as to avoid double-naming models from a different script

EBMod10 <- glmer(cbind(dead,live) ~ temperature*treatment*day + (1|strain), family = binomial(link="logit"), data = EBData)
summary(EBMod10)

#and model checks:
plot(fitted(EBMod10), residuals(EBMod10))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod10))
qqline(residuals(EBMod10), col = "red")

# Deviance and AIC
deviance(EBMod10)
#1744.379
AIC(EBMod10)
#2914.154

# Overdispersion check
overdispersion <- deviance(EBMod10) / df.residual(EBMod10)
overdispersion
#hmm... could be too overdispersed...


###################

##### Now for the model with strain as a random effect and ID as random effect:


EBMod11 <- glmer(cbind(dead,live) ~ temperature*treatment*day + (1|strain) + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod11)

#and model checks:
plot(fitted(EBMod11), residuals(EBMod11))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod11))
qqline(residuals(EBMod11), col = "red")

# Deviance and AIC
deviance(EBMod11)
#1306.417
AIC(EBMod11)
#2638.688

# Overdispersion check
overdispersion <- deviance(EBMod11) / df.residual(EBMod11)
overdispersion
#hmm... a bit better.
#5.680073



#########################################

#### Trying with scaled variables:


### Well fo course, I need to scale the variables:

### Code to scale variables ###
EBData <- EBData %>%
  mutate(
    temperature_scaled = scale(temperature),
    day_scaled = scale(day)
  )



###### Model with only strain as random effect ###########

EBMod12 <- glmer(cbind(dead,live) ~ temperature_scaled*treatment*day_scaled + (1|strain), family = binomial(link="logit"), data = EBData)
summary(EBMod12)


#and model checks:
plot(fitted(EBMod12), residuals(EBMod12))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod12))
qqline(residuals(EBMod12), col = "red")

# Deviance and AIC
deviance(EBMod12)
#1744.379
AIC(EBMod12)
#2914.154

# Overdispersion check
overdispersion <- deviance(EBMod12) / df.residual(EBMod12)
overdispersion
#hmm... a bit better.
#7.551424



##### Model with strain and ID as random effect ############

EBMod13 <- glmer(cbind(dead,live) ~ temperature_scaled*treatment*day_scaled + (1|strain) + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod13)

#and model checks:
plot(fitted(EBMod13), residuals(EBMod13))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod13))
qqline(residuals(EBMod13), col = "red")

# Deviance and AIC
deviance(EBMod13)
#1306.418
AIC(EBMod13)
#2638.688

# Overdispersion check
overdispersion <- deviance(EBMod13) / df.residual(EBMod13)
overdispersion
#5.680079



#### Calculating proportion predictions based on this model:



### First I need to find a way to understand what the unscaled values are, and which
#scaled values I am interested in:



### 
# day 7 scaled = 
(50 - day_mean)/day_sd
#day 0 = -2.08584
#day 7 = -1.599144
#day 14 = -1.112448
#day 21 = -0.625752
#day 28 = -0.139056
#day 35 = 0.34764
#day 42 = 0.834336
#day 49 = 1.321032
#day 50 = 1.39056
#day 56 = 1.807728


#temperature 1C scaled = 
(4-temperature_mean)/temperature_sd
#1C = -0.9979145
#4C = 0.9979145



### averaging across strain and Replicate ID:

calculate_prediction_average <- function(temperature_scaled, treatment, day_scaled) {
  # Fixed effects from your model output
  intercept <- -2.41677
  beta_temp <- 0.07847
  beta_treat <- 1.05444
  beta_day <- 0.02507
  interaction_temp_treat <- 0.08209
  interaction_temp_day <- 0.05541
  interaction_treat_day <- 0.49444
  interaction_all <- -0.03663
  
  # Calculate the linear predictor (eta) without random effects
  eta <- intercept +
    beta_temp * temperature_scaled +
    beta_treat * (treatment == "Dark") +
    beta_day * day_scaled +
    interaction_temp_treat * temperature_scaled * (treatment == "Dark") +
    interaction_temp_day * temperature_scaled * day_scaled +
    interaction_treat_day * (treatment == "Dark") * day_scaled +
    interaction_all * temperature_scaled * (treatment == "Dark") * day_scaled
  
  # Convert linear predictor to probability
  p <- 1 / (1 + exp(-eta))
  return(p)
}

#day 0 = -2.08584
#day 7 = -1.599144
#day 14 = -1.112448
#day 21 = -0.625752
#day 28 = -0.139056
#day 35 = 0.34764
#day 42 = 0.834336
#day 49 = 1.321032
#day 50 = 1.39056
#day 56 = 1.807728

# Finding average proportions of dead cells:
#1C:
#Dark:
calculate_prediction_average(-0.9979145, "Dark", 1.39056)

#Light:
calculate_prediction_average(-0.9979145, "Light", 1.39056)


#4C:
#Dark:
calculate_prediction_average(0.9979145, "Dark", 1.39056)

#Light:
calculate_prediction_average(0.9979145, "Light", 1.39056)




################################################################################

#Finally, a post-hoc model to investigate the differences between strains:


EBMod5b <- glmer(cbind(dead,live) ~ temperature*treatment*day + strain + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod5b)
# comes up with warnings again


plot(fitted(EBMod5b), residuals(EBMod5b))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod5b))
qqline(residuals(EBMod5b), col = "red")

# Deviance and AIC
deviance(EBMod5b)
#1307.928
AIC(EBMod5b)
#2625.607

# Overdispersion check
overdispersion <- deviance(EBMod5b) / df.residual(EBMod5b)
overdispersion
#5.761799


###############################

# Let's try with the scaled variables:

EBMod5c <- glmer(cbind(dead,live) ~ temperature_scaled*treatment*day_scaled + strain + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod5c)
#no warnings


plot(fitted(EBMod5c), residuals(EBMod5c))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod5c))
qqline(residuals(EBMod5c), col = "red")

# Deviance and AIC
deviance(EBMod5c)
#1307.927
AIC(EBMod5c)
#2625.607

# Overdispersion check
overdispersion <- deviance(EBMod5c) / df.residual(EBMod5c)
overdispersion
#5.761794








