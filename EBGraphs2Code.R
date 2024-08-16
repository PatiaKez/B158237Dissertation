#### Graphs for EB Attempt n2 #####

library(ggplot2)
library(dplyr)
library(gridExtra)



# Fit regression models for each combination of strain, treatment, and temperature
regression_models <- EBData %>%
  group_by(strain, treatment, temperature) %>%
  do(model = lm(prop_dead ~ day, data = .))

# Create predictions for plotting trend lines
predictions <- regression_models %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pred_data = list(
    data.frame(day = seq(min(EBData$day), max(EBData$day), length.out = 100),
               prop_dead = predict(model, newdata = data.frame(day = seq(min(EBData$day), max(EBData$day), length.out = 100)))
    )
  )) %>%
  unnest(pred_data)

# Calculate mean and standard error for plotting points
summary_data <- EBData %>%
  group_by(strain, treatment, temperature, day) %>%
  summarize(
    mean_prop_dead = mean(prop_dead, na.rm = TRUE),
    sd_prop_dead = sd(prop_dead, na.rm = TRUE),
    n = n(),
    se_prop_dead = sd_prop_dead / sqrt(n),
    .groups = 'drop'
  )

# Filter and plot for temperature 0.5
summary_temp_05 <- filter(summary_data, temperature == 0.5)
predictions_temp_05 <- filter(predictions, temperature == 0.5)

plot_temp_05 <- ggplot() +
  geom_line(data = predictions_temp_05, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_05, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_05, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    title = "Mean Proportion of Dead Cells Over Time at 1C Temperature",
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
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
  geom_line(data = predictions_temp_4, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_4, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_4, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    title = "Mean Proportion of Dead Cells Over Time at 4C Temperature",
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
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
grid.arrange(plot_temp_05, plot_temp_4, ncol = 2)


###########################


# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Filter and plot for temperature 0.5
summary_temp_05 <- filter(summary_data, temperature == 1)
predictions_temp_05 <- filter(predictions, temperature == 1)

plot_temp_05 <- ggplot() +
  geom_line(data = predictions_temp_05, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_05, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_05, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Day",
    y = "Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# Filter and plot for temperature 4
summary_temp_4 <- filter(summary_data, temperature == 4)
predictions_temp_4 <- filter(predictions, temperature == 4)

plot_temp_4 <- ggplot() +
  geom_line(data = predictions_temp_4, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_4, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_4, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Day",
    y = "Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95") 
  )

# Combine plots into one display
grid.arrange(plot_temp_05, plot_temp_4, ncol = 2)

