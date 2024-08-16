######### Analysing only the Dark treatment data #############

# Load the same full dataset

#loading in the data, that might be useful
#EBData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/MScDiatomRCode/EBData.csv", header = TRUE)

#And create the proportions column

## First I will create a new column and transform my percentage into proportion dead:
#EBData$prop_dead <- EBData$perc_dead/100

#####

# Now I want to subset the full dataframe to only include the dark treatment

EBData_dark <- EBData %>% filter(treatment == "Dark")

head(EBData_dark)

#lovely, that seems to have worked

#######

## Now I need to repeat my stats on this data

#First I will visualise the data:

ggplot(EBData_dark, aes(x = day, y = prop_dead, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Proportion of Dead Cells Over Time",
       x = "Timepoint (weeks)",
       y = "Proportion of Dead Cells") +
  theme_minimal()

### Well, in retrospect, this was a little pointless... I have already visualised the data, 
#i might have less graphs, but not like the data will change.

#ok, so let's look at the histogram of the data distribution:

hist(EBData_dark$prop_dead)
#hmm... that appears to be approaching normality but still with a left-skew 

#I will try repeat the binomial model I had in the first part of analysis:
EBDMod1 <- glm(cbind(dead,live) ~ temperature*day + strain, family = binomial(link="logit"), data = EBData_dark)
summary(EBDMod1)

plot(EBDMod1)


plot(fitted(EBDMod1), residuals(EBDMod1))
abline(h = 0, col = "red")

qqnorm(residuals(EBDMod1))
qqline(residuals(EBDMod1), col = "red")

# Deviance and AIC
deviance(EBDMod1)
#918.0111
AIC(EBDMod1)
#1526.944

# Overdispersion check
overdispersion <- deviance(EBDMod1) / df.residual(EBDMod1)
overdispersion
#hmm... could be too overdispersed...
#8.196528


######################

EBDMod2 <- glmer(cbind(dead,live) ~ temperature*day + (1|strain), family = binomial(link="logit"), data = EBData_dark)
summary(EBDMod2)
#comes up with warnings

plot(fitted(EBDMod2), residuals(EBDMod2))
abline(h = 0, col = "red")

qqnorm(residuals(EBDMod2))
qqline(residuals(EBDMod2), col = "red")

# Deviance and AIC
deviance(EBDMod2)
#918.0383
AIC(EBDMod2)
#1551.719

# Overdispersion check
overdispersion <- deviance(EBDMod2) / df.residual(EBDMod2)
overdispersion
#hmm... could be too overdispersed...
#7.982942


############################################

EBDMod3 <- glmer(cbind(dead,live) ~ temperature*day + (1|strain) + (1|ID), family = binomial(link="logit"), data = EBData_dark)
summary(EBDMod3)
## comes up with scaling warnings

plot(fitted(EBDMod3), residuals(EBDMod3))
abline(h = 0, col = "red")

qqnorm(residuals(EBDMod3))
qqline(residuals(EBDMod3), col = "red")

# Deviance and AIC
deviance(EBDMod3)
#749.3244
AIC(EBDMod3)
#1455.187

# Overdispersion check
overdispersion <- deviance(EBDMod3) / df.residual(EBDMod3)
overdispersion
#hmm... could be too overdispersed...
#6.573021



############################################

#### Trying with scaled variables:


### Well fo course, I need to scale the variables:

### Code to scale variables ###
EBData_dark <- EBData_dark %>%
  mutate(
    temperature_scaled = scale(temperature),
    day_scaled = scale(day)
  )


#########

EBDMod4 <- glmer(cbind(dead,live) ~ temperature_scaled*day_scaled + (1|strain), family = binomial(link="logit"), data = EBData_dark)
summary(EBDMod4)
#no warnings


plot(fitted(EBDMod4), residuals(EBDMod4))
abline(h = 0, col = "red")

qqnorm(residuals(EBDMod4))
qqline(residuals(EBDMod4), col = "red")

# Deviance and AIC
deviance(EBDMod4)
#918.0383
AIC(EBDMod4)
#1551.719

# Overdispersion check
overdispersion <- deviance(EBDMod4) / df.residual(EBDMod4)
overdispersion
#hmm... could be too overdispersed...
#7.982942



################################################

EBDMod5 <- glmer(cbind(dead,live) ~ temperature_scaled*day_scaled + (1|strain) + (1|ID), family = binomial(link="logit"), data = EBData_dark)
summary(EBDMod5)
#no warnings

plot(fitted(EBDMod5), residuals(EBDMod5))
abline(h = 0, col = "red")

qqnorm(residuals(EBDMod5))
qqline(residuals(EBDMod5), col = "red")

# Deviance and AIC
deviance(EBDMod5)
#749.3245
AIC(EBDMod5)
#1455.187

# Overdispersion check
overdispersion <- deviance(EBDMod5) / df.residual(EBDMod5)
overdispersion
#hmm... could be too overdispersed...
#6.573022



###### Ok, so it looks like the models with both strain and replicate ID as random effects
## give the best fit (slightly worse QQ Plot, bot better dispersion and AIC)/
## Both give similar results, but the unscaled model gives a very large eigenvalue warning
## and suggests to re-scale variables.
## I need to choose between these two.



#####################################################


#### So why don't I try to make some plots of this

#What I think I should do is get means of each replicate trio and the standard errors

#so first I need to group the data by light, temperature and strain:

#load the necessary libraries:

library(dplyr)
library(tidyr)


# Group the data by strain, temperature, and timepoint, then calculate mean and standard error

EBDGrouped <- EBData_dark %>%
  group_by(strain, temperature, day) %>%
  summarise(
    mean_proportion_dead = mean(prop_dead, na.rm = TRUE),
    se_proportion_dead = sd(prop_dead, na.rm = TRUE) / sqrt(n())
  )


############################################################


# Fit regression models for each combination of strain, treatment, and temperature
regression_models <- EBData_dark %>%
  group_by(strain, temperature) %>%
  do(model = lm(prop_dead ~ day, data = .))

# Create predictions for plotting trend lines
predictions <- regression_models %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pred_data = list(
    data.frame(day = seq(min(EBData_dark$day), max(EBData_dark$day), length.out = 100),
               prop_dead = predict(model, newdata = data.frame(day = seq(min(EBData_dark$day), max(EBData_dark$day), length.out = 100)))
    )
  )) %>%
  unnest(pred_data)

# Calculate mean and standard error for plotting points
summary_data <- EBData_dark %>%
  group_by(strain, treatment, temperature, day) %>%
  summarize(
    mean_prop_dead = mean(prop_dead, na.rm = TRUE),
    sd_prop_dead = sd(prop_dead, na.rm = TRUE),
    n = n(),
    se_prop_dead = sd_prop_dead / sqrt(n),
    .groups = 'drop'
  )

# Filter and plot for temperature 1
summary_temp_1 <- filter(summary_data, temperature == 1)
predictions_temp_1 <- filter(predictions, temperature == 1)

plot_temp_1 <- ggplot() +
  geom_line(data = predictions_temp_1, aes(x = day, y = prop_dead, color = strain), size = 1.2) +
  geom_point(data = summary_temp_1, aes(x = day, y = mean_prop_dead, color = strain), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_1, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    title = "Mean Proportion of Dead Cells Over Time at 1C Temperature",
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# Filter and plot for temperature 4
summary_temp_4 <- filter(summary_data, temperature == 4)
predictions_temp_4 <- filter(predictions, temperature == 4)

plot_temp_4 <- ggplot() +
  geom_line(data = predictions_temp_4, aes(x = day, y = prop_dead, color = strain), size = 1.2) +
  geom_point(data = summary_temp_4, aes(x = day, y = mean_prop_dead, color = strain,), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_4, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    title = "Mean Proportion of Dead Cells Over Time at 4C Temperature",
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# Combine plots into one display
grid.arrange(plot_temp_1, plot_temp_4, ncol = 2)


###########################


# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




