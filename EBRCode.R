### MSc Antarctic Diatom Project - R Code ###

#where to even start...
##I will load my necessary libraries...

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)

#loading in the data, that might be useful
EBData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/MScDiatomRCode/EBData.csv", header = TRUE)

names(EBData)[1] <- "ID"

## First I will create a new column and transform my percentage into proportion dead:

EBData$prop_dead <- EBData$perc_dead/100

#fist plot I can think of...
plot(perc_dead~day, data =EBData)

#ok... not as terrible as I thought, I now need to find a way to separate it by strain...
#and by treatment


## I want to visualise the data  bit:

# Create the plot
ggplot(EBData, aes(x = day, y = prop_dead, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Proportion of Dead Cells Over Time",
       x = "Timepoint (weeks)",
       y = "Proportion of Dead Cells") +
  theme_minimal()

#so it looks like there might be an effect of the darkness on proportion of dead cells

#Let's try another one
ggplot(EBData, aes(x = as.factor(timepoint), y = prop_dead, fill = strain)) +
  geom_boxplot() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Proportion of Dead Cells at Different Timepoints",
       x = "Timepoint (weeks)",
       y = "Proportion of Dead Cells") +
  theme_minimal()

#hmm... seems like theres quite a lot of variance in each measurement


## Now onto statistical models

install.packages("lmerTest")
library(lmerTest)

EBMod1 <- lm(perc_dead ~ strain + temperature + treatment + day, data = EBData)
summary(EBMod1)


#let's see model diagnostics?
plot(EBMod1)

#not too bad for a preliminary model actually

#ok, that was my preliminary test model...
#now on to more complete, sophisticated models with random effects and interactions

EBMod2 <- lm(perc_dead ~ strain + temperature*treatment*day, data = EBData)
summary(EBMod2)
#ok...

#model diagnostics:
plot(EBMod2)
#hmmm... QQ plot is not so good anymore

#Let's try another one:

EBMod3 <- lmer(perc_dead ~ temperature*treatment*day + (1|strain), data = EBData)
summary(EBMod3)

#and model diagnostics:
plot(EBMod3)

#ahh, wait, for mixed models I need a different package to get diagnostic plots...
#which one again?

EBMod4 <- lm(perc_dead ~ strain*temperature*treatment*day, data = EBData)
summary(EBMod4)
#yeah... I think that's entirely too much

plot(EBMod4)

#BBFAMCmod <- glmer(beefav ~  mc  + (1|site) + (1|surveyround), data = flowersclass, family = poisson)
#summary(BBFAMCmod)

#####################


### Maybe next I should try to make actual graphs of the effects, so that I can interpret and 
#improve the stats. i.e. choose which models are best representative of the data.

hist(EBData$perc_dead)
# that looks like it might be a poisson

mean(EBData$perc_dead)

var(EBData$perc_dead)

# for a poisson the mean = variance

##########################

install.packages("emmeans")
library(emmeans)



## Let's step back...
#I have percentage data as my response variable, this is basically equivalent to 
#proportional data.

## I can analyse proportional data with a binomial distribution and a logit link
#this will bound my analysis at 0 and 1.



## Next I will fit a binomial regression with a logit link:

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


# Cook's distance
cooksd <- cooks.distance(EBMod5)
plot(cooksd, type = "h")
abline(h = 4 / nrow(EBData), col = "red") # Rule of thumb for identifying influential points
#??


# Install and load the car package if not already installed
install.packages("car")
library(car)

# Calculate VIF
vif(EBMod5, type = "predictor")




## Trying to add a replicate-level random effect

EBMod5b <- glmer(cbind(dead,live) ~ temperature*treatment*day + strain + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod5b)

##once again getting warnings... brilliant.

plot(fitted(EBMod5b), residuals(EBMod5b))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod5b))
qqline(residuals(EBMod5b), col = "red")

# Deviance and AIC
deviance(EBMod5b)
AIC(EBMod5b)

# Overdispersion check
overdispersion <- deviance(EBMod5b) / df.residual(EBMod5b)
overdispersion
#hmm... could be too overdispersed...


# Cook's distance
cooksd <- cooks.distance(EBMod5b)
plot(cooksd, type = "h")
abline(h = 4 / nrow(EBData), col = "red") # Rule of thumb for identifying influential points
#??

### hmmm, perhaps a bit better on overdispersion, but otherwise not that different
# but let's check:

AO56 <- anova(EBMod5,EBMod5b)
print(AO56)

anova(EBMod5, EBMod5b, test = "Chisq")

################################

### Code to scale variables ###
EBData <- EBData %>%
  mutate(
    temperature_scaled = scale(temperature),
    day_scaled = scale(day)
  )

#####

EBData$day <- as.numeric(EBData$day)

control <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

EBMod6 <- glmer(cbind(dead,live) ~ temperature*treatment*day +(1|strain), family = binomial(link="logit"), data = EBData, control = control)
summary(EBMod6)

library(car)
vif(EBMod6)


EBMod6b <- glmer(cbind(dead,live) ~ temperature*treatment*day +(1|strain) + (1|ID), family = binomial(link="logit"), data = EBData, control = control)
summary(EBMod6)



####

EBMod6a <- glmer(cbind(dead,live) ~ temperature_scaled*treatment*day_scaled +(1|strain) + (1|ID), family = binomial(link="logit"), data = EBData)
summary(EBMod6a)

plot(EBMod6a)

plot(fitted(EBMod6a), residuals(EBMod6a))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod6a))
qqline(residuals(EBMod6a), col = "red")

# Deviance and AIC
deviance(EBMod6a)
AIC(EBMod6a)

# Overdispersion check
overdispersion <- deviance(EBMod6a) / df.residual(EBMod6a)
overdispersion
#hmm... could be too overdispersed...


# Cook's distance
cooksd <- cooks.distance(EBMod5)
plot(cooksd, type = "h")
abline(h = 4 / nrow(EBData), col = "red") # Rule of thumb for identifying influential points
#??


###### EBMod5 Post-hoc Tukey test #########

##########################
library(multcomp)

# Subset data for different combinations of temperature and light
df_temp4_light <- EBData %>% filter(temperature == 4 & treatment == "Control")
df_temp4_dark <- EBData %>% filter(temperature == 4 & treatment == "Dark")
df_temp1_light <- EBData %>% filter(temperature == 1 & treatment == "Control")
df_temp1_dark <- EBData %>% filter(temperature == 1 & treatment == "Dark")

# Fit models for each subset
model_temp4_light <- glm(cbind(dead,live) ~ day + strain, family = binomial(link = "logit"), data = df_temp4_light)
model_temp4_dark <- glm(cbind(dead,live) ~ day + strain, family = binomial(link = "logit"), data = df_temp4_dark)
model_temp1_light <- glm(cbind(dead,live) ~ day + strain, family = binomial(link = "logit"), data = df_temp1_light)
model_temp1_dark <- glm(cbind(dead,live) ~ day + strain, family = binomial(link = "logit"), data = df_temp1_dark)

summary(model_temp4_light)
summary(model_temp4_dark)
summary(model_temp1_light)
summary(model_temp1_dark)

##########################


# Compute estimated marginal means
emmeans_model <- emmeans(EBMod5, ~ day * temperature * treatment | strain)

# Perform pairwise comparisons with Tukey adjustment
pairwise_comparisons <- contrast(emmeans_model, method = "pairwise", adjust = "tukey")
summary(pairwise_comparisons)


library(ggplot2)

# Example visualization of interaction effects
#EBData$day <- as.factor(EBData$day)  # Convert day to factor if needed
ggplot(EBData, aes(x = day, y = prop_dead, color = as.factor(temperature), linetype = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ strain) +
  labs(title = "Interaction Effects", x = "Day", y = "Proportion of Dead Cells") +
  theme_minimal()


##########################


####### Post Hoc Tukey test for Model EBMod6 ##################

# Compute estimated marginal means
emmeans_model6 <- emmeans(EBMod6, ~ day * temperature * treatment)

# Perform pairwise comparisons with Tukey adjustment
pairwise_comparisons6 <- contrast(emmeans_model6, method = "pairwise", adjust = "tukey")
summary(pairwise_comparisons)

##
library(ggplot2)

ggplot(EBData, aes(x = day, y = (dead / (dead + live)), color = as.factor(temperature), linetype = treatment)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ strain) +
  labs(title = "Interaction Effects", x = "Day", y = "Proportion of Dead Cells") +
  theme_minimal()





############################################

#### Trying some further models #####


EBMod7 <- glm(cbind(dead,live) ~ temperature*treatment*day*strain, family = binomial(link="logit"), data = EBData)
summary(EBMod7)

plot(EBMod7)
#hmm... nightmare to interpret though...

##### Quasi-binomial distribution:

EBMod8 <- glm(cbind(dead, live) ~ temperature * treatment * day + strain, 
              family = quasibinomial(link = "logit"), data = EBData)
summary(EBMod8)


# Residual plots
par(mfrow = c(2, 2))
plot(EBMod6)

# Residuals vs fitted values
plot(fitted(EBMod7), residuals(EBMod7))
abline(h = 0, col = "red")

# Q-Q plot
qqnorm(residuals(EBMod7))
qqline(residuals(EBMod7), col = "red")

# Deviance and AIC
deviance(EBMod7)
AIC(EBMod7)

# Cook's distance
cooksd <- cooks.distance(EBMod7)
plot(cooksd, type = "h")
abline(h = 4 / nrow(EBData), col = "red")

# Overdispersion check for quasi-binomial model
overdispersion <- deviance(EBMod7) / df.residual(EBMod7)
overdispersion
# well... overdispersion is not reduced, might not be a solution...

