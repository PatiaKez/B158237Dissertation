###### Latitude graph

library(ggplot2)
library(dplyr)
library(tidyverse)

#loading in the data, that might be useful
LatData <- read.csv("C:/Users/User/Documents/Edi Masters/Master's Research Project/MScDiatomRCode/MortalityByLatitudeData.csv", header = TRUE)

#
ggplot(LatData, aes(x = latitude, y = perc_dead, color = as.factor(Temperature), group = Temperature)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("blue", "red"), labels = c("1°C", "4°C")) +
  labs(
    title = "Percentage of Dead Cells at Different Latitudes",
    x = "Latitude",
    y = "Percentage of Dead Cells",
    color = "Temperature"
  ) +
  theme_minimal()



## ok not bad, now to make it better:


ggplot(LatData, aes(x = latitude, y = perc_dead, color = as.factor(Temperature), group = Temperature)) +
  geom_line(size = 1) +  
  geom_point(size = 2.5) +  
  scale_color_manual(values = c("blue", "red"), labels = c("1°C", "4°C")) +  
  scale_x_continuous(
    breaks = 1:8,  
    labels = c("66 (0 days)", "67.4 (7 days)", "68.2 (14 days)", "69.1 (21 days)", "70 (28 days)", "71 (35 days)", "72.1 (42 days)", "73.2 (49 days)")  # Custom labels
  ) +
  coord_cartesian(ylim = c(0, 50)) +  
  labs(
    x = "Latitude and Corresponding Maximum Continuous Dark Duration",
    y = "Mortality (% dead cells)",
    color = "Temperature"
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13), 
    axis.text.y = element_text(size = 13),  # Enlarge y-axis labels
    axis.title.x = element_text(size = 15),  # Enlarge x-axis title
    axis.title.y = element_text(size = 15)   # Enlarge y-axis title
  )





