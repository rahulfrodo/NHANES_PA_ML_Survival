setwd("D:/Final_Year_Project")
rm(list = ls())    
Sys.setenv(LANGUAGE = "en")
gc() # free unused memory

library(ggplot2)
library(scales)
library(patchwork)

# Load data
c_index_df <- read.csv("c_index_finalPlot.csv")
c_index_df$Split <- 1:nrow(c_index_df)

###################Cox
# Compute means
mean_c_index_cox_MeanAct <- mean(c_index_df$c_index_cox_MeanAct)
mean_c_index_cox_FPC <- mean(c_index_df$c_index_cox_FPC)


###
sd(c_index_df$c_index_cox_MeanAct)
sd(c_index_df$c_index_Coxtime_FPC)

sd(c_index_df$c_index_rsf_MeanAct)
sd(c_index_df$c_index_rsf_FPC)

sd(c_index_df$c_index_gbsm_MeanAct)
sd(c_index_df$c_index_gbsm_FPC)

sd(c_index_df$c_index_deepsurv_MeanAct)
sd(c_index_df$c_index_deepsurv_FPC)


###


# Prepare colors
default_colors <- hue_pal()(2)
model_colors <- c(
  "Cox Model(MeanAct)" = default_colors[1],
  "Cox Model(FPC)" = default_colors[2]
)

# Main plot
main_plot <- ggplot() +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_cox_MeanAct, 
                                   color = "Cox Model(MeanAct)")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_cox_FPC, 
                                   color = "Cox Model(FPC)")) +
  
  geom_point(data = c_index_df, aes(x = Split, y = c_index_cox_MeanAct, 
                                    color = "Cox Model(MeanAct)")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_cox_FPC, 
                                    color = "Cox Model(FPC)")) +
  
  geom_hline(aes(yintercept = mean_c_index_cox_MeanAct, 
                 color = "Cox Model(MeanAct)"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept =mean_c_index_cox_FPC, 
                 color = "Cox Model(FPC)"), linetype = "dashed", size = 1) +
  
  scale_color_manual(values = model_colors) +
  scale_y_continuous(
    breaks = seq(0.725, 0.85, by = 0.0125),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  coord_cartesian(ylim = c(0.725, 0.85)) +
  
  # add space on both ends of x-axis
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Split Number",
    y = "C-index",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Bottom caption with updated means
mean_label <- paste(
  "Mean C-index | Cox Model(MeanAct):", round(mean_c_index_cox_MeanAct, 5),
  "| Cox Model(FPC):", round(mean_c_index_cox_FPC, 5)
)

bottom_caption <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = mean_label, size = 4, 
           hjust = 0.5) +
  theme_void()

# Combine all
final_plot_cox <- main_plot / bottom_caption + 
  plot_layout(heights = c(1, 0.1))
print(final_plot_cox)

###################RSF
# Compute means
mean_c_index_rsf_MeanAct <- mean(c_index_df$c_index_rsf_MeanAct)
mean_c_index_rsf_FPC <- mean(c_index_df$c_index_rsf_FPC)


# Prepare colors
default_colors <- hue_pal()(2)
model_colors <- c(
  "RSF Model(MeanAct)" = default_colors[1],
  "RSF Model(FPC)" = default_colors[2]
)

# Main plot
main_plot <- ggplot() +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_rsf_MeanAct, 
                                   color = "RSF Model(MeanAct)")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_rsf_FPC, 
                                   color = "RSF Model(FPC)")) +
  
  geom_point(data = c_index_df, aes(x = Split, y = c_index_rsf_MeanAct, 
                                    color = "RSF Model(MeanAct)")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_rsf_FPC, 
                                    color = "RSF Model(FPC)")) +
  
  geom_hline(aes(yintercept = mean_c_index_rsf_MeanAct, 
                 color = "RSF Model(MeanAct)"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept =mean_c_index_rsf_FPC, 
                 color = "RSF Model(FPC)"), linetype = "dashed", size = 1) +
  
  scale_color_manual(values = model_colors) +
  scale_y_continuous(
    breaks = seq(0.725, 0.85, by = 0.0125),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  coord_cartesian(ylim = c(0.725, 0.85)) +
  
  # add space on both ends of x-axis
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Split Number",
    y = "C-index",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Bottom caption with updated means
mean_label <- paste(
  "Mean C-index | RSF Model(MeanAct):", round(mean_c_index_rsf_MeanAct, 5),
  "| RSF Model(FPC):", round(mean_c_index_rsf_FPC, 5)
)

bottom_caption <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = mean_label, size = 4, 
           hjust = 0.5) +
  theme_void()

# Combine all
final_plot_rsf <- main_plot / bottom_caption + 
  plot_layout(heights = c(1, 0.1))
print(final_plot_rsf)

####################GBSM
# Compute means
mean_c_index_gbsm_MeanAct <- mean(c_index_df$c_index_gbsm_MeanAct)
mean_c_index_gbsm_FPC <- mean(c_index_df$c_index_gbsm_FPC)


# Prepare colors
default_colors <- hue_pal()(2)
model_colors <- c(
  "GBSM Model(MeanAct)" = default_colors[1],
  "GBSM Model(FPC)" = default_colors[2]
)

# Main plot
main_plot <- ggplot() +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_gbsm_MeanAct, 
                                   color = "GBSM Model(MeanAct)")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_gbsm_FPC, 
                                   color = "GBSM Model(FPC)")) +
  
  geom_point(data = c_index_df, aes(x = Split, y = c_index_gbsm_MeanAct, 
                                    color = "GBSM Model(MeanAct)")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_gbsm_FPC, 
                                    color = "GBSM Model(FPC)")) +
  
  geom_hline(aes(yintercept = mean_c_index_gbsm_MeanAct, 
                 color = "GBSM Model(MeanAct)"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept =mean_c_index_gbsm_FPC, 
                 color = "GBSM Model(FPC)"), linetype = "dashed", size = 1) +
  
  scale_color_manual(values = model_colors) +
  scale_y_continuous(
    breaks = seq(0.725, 0.85, by = 0.0125),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  coord_cartesian(ylim = c(0.725, 0.85)) +
  
  # add space on both ends of x-axis
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Split Number",
    y = "C-index",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Bottom caption with updated means
mean_label <- paste(
  "Mean C-index | GBSM Model(MeanAct):", round(mean_c_index_gbsm_MeanAct, 5),
  "| GBSM Model(FPC):", round(mean_c_index_gbsm_FPC, 5)
)

bottom_caption <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = mean_label, size = 4, 
           hjust = 0.5) +
  theme_void()

# Combine all
final_plot_gbsm <- main_plot / bottom_caption + 
  plot_layout(heights = c(1, 0.1))
print(final_plot_gbsm)

####################Deepsurv
# Compute means
mean_c_index_deepsurv_MeanAct <- mean(c_index_df$c_index_deepsurv_MeanAct)
mean_c_index_deepsurv_FPC <- mean(c_index_df$c_index_deepsurv_FPC)


# Prepare colors
default_colors <- hue_pal()(2)
model_colors <- c(
  "Deepsurv Model(MeanAct)" = default_colors[1],
  "Deepsurv Model(FPC)" = default_colors[2]
)

# Main plot
main_plot <- ggplot() +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_deepsurv_MeanAct, 
                                   color = "Deepsurv Model(MeanAct)")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_deepsurv_FPC, 
                                   color = "Deepsurv Model(FPC)")) +
  
  geom_point(data = c_index_df, aes(x = Split, y = c_index_deepsurv_MeanAct, 
                                    color = "Deepsurv Model(MeanAct)")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_deepsurv_FPC, 
                                    color = "Deepsurv Model(FPC)")) +
  
  geom_hline(aes(yintercept = mean_c_index_deepsurv_MeanAct, 
                 color = "Deepsurv Model(MeanAct)"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept =mean_c_index_deepsurv_FPC, 
                 color = "Deepsurv Model(FPC)"), linetype = "dashed", size = 1) +
  
  scale_color_manual(values = model_colors) +
  scale_y_continuous(
    breaks = seq(0.725, 0.85, by = 0.0125),
    labels = scales::label_number(accuracy = 0.001)
  ) +
  coord_cartesian(ylim = c(0.725, 0.85)) +
  
  # add space on both ends of x-axis
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Split Number",
    y = "C-index",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Bottom caption with updated means
mean_label <- paste(
  "Mean C-index | Deepsurv Model(MeanAct):", round(mean_c_index_deepsurv_MeanAct, 5),
  "| Deepsurv Model(FPC):", round(mean_c_index_deepsurv_FPC, 5)
)

bottom_caption <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = mean_label, size = 4, 
           hjust = 0.5) +
  theme_void()

# Combine all
final_plot_deepsurv <- main_plot / bottom_caption + 
  plot_layout(heights = c(1, 0.1))
print(final_plot_deepsurv)
dev.off()

###########################Combine all plots into a panel


library(patchwork)

# Arrange in 2 columns and 3 rows
panel_plot <- (
  final_plot_cox | final_plot_rsf
) / (
  final_plot_gbsm | final_plot_deepsurv
) 

# Add a global title and optional tags
panel_plot <- panel_plot + 
  plot_annotation(
    title = "Comparison of Model Performance Using C-index Across 100 Splits",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ),
  )

# Print the complete plot
print(panel_plot)
