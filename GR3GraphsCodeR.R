######### Graphs for Growth Rate #######################

## Load necessary libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Calculate growth rates (example calculation, adjust based on your actual data)
GRData1 <- GRData %>%
  group_by(strain, treatment, temperature, replicate) %>%
  arrange(day) %>%
  mutate(growth_rate = growthrateB) %>%
  filter(!is.na(growth_rate))

# Fit the linear models to the growth rates
regression_models <- GRData1 %>%
  group_by(strain, treatment, temperature) %>%
  do(model = lm(growth_rate ~ day, data = .))

# Create predictions for plotting trend lines
predictions <- regression_models %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pred_data = list(
    data.frame(day = seq(min(GRData1$day), max(GRData1$day), length.out = 100),
               growth_rate = predict(model, newdata = data.frame(day = seq(min(GRData1$day), max(GRData1$day), length.out = 100)))
    )
  )) %>%
  unnest(pred_data)

# Calculate mean and standard error for plotting points
summary_data <- GRData1 %>%
  group_by(strain, treatment, temperature, day) %>%
  summarize(
    mean_growth_rate = mean(growth_rate, na.rm = TRUE),
    sd_growth_rate = sd(growth_rate, na.rm = TRUE),
    n = n(),
    se_growth_rate = sd_growth_rate / sqrt(n),
    .groups = 'drop'
  )

# Function to create plots by temperature
create_plot <- function(summary_data, predictions, temperature, title) {
  summary_filtered <- filter(summary_data, temperature == !!temperature)
  predictions_filtered <- filter(predictions, temperature == !!temperature)
  
  ggplot() +
    geom_line(data = predictions_filtered, aes(x = day, y = growth_rate, color = strain, linetype = treatment), size = 1.2) +
    geom_point(data = summary_filtered, aes(x = day, y = mean_growth_rate, color = strain, shape = treatment), size = 4, stroke = 1.5) +
    geom_errorbar(data = summary_filtered, aes(x = day, ymin = mean_growth_rate - se_growth_rate, ymax = mean_growth_rate + se_growth_rate), width = 0.3) +
    labs(
      x = "Day",
      y = "Growth Rate",
      color = "Strain",
      shape = "Light Treatment",
      linetype = "Light Treatment"
    )+ scale_color_manual(values = colorblind_palette) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
}

# Create plots for each temperature
plot_temp_05 <- create_plot(summary_data, predictions, 0.5, "Growth Rate at 1C Temperature")
plot_temp_4 <- create_plot(summary_data, predictions, 4, "Growth Rate at 4C Temperature")

# Combine plots into one display
grid.arrange(plot_temp_05, plot_temp_4, ncol = 2, nrow = 1)

