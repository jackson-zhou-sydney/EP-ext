
# Analysis of results

source("setup.R")

l1_cdf <- data.frame(method = character(),
                     init = double(),
                     j = double(),
                     l1 = double())

other_cdf <- data.frame(method = character(),
                        init = double(),
                        time = double(),
                        mem = double())

res_files <- list.files("results/")

for (res_file in res_files) {
  if (res_file != "mcmc.RData") {
    load(paste0("results/", res_file))
    l1_cdf <- rbind(l1_cdf, l1_df)
    other_cdf <- rbind(other_cdf, other_df)
  }
}

l1_cdf_summary <- l1_cdf |> group_by(method, init) |>  
  dplyr::summarise(mean_l1 = mean(l1))

mean_cdf <- cbind(other_cdf, mean_l1 = l1_cdf_summary$mean_l1) |> 
  group_by(method) |> 
  dplyr::summarise(mean_l1 = mean(mean_l1),
                   time = mean(time),
                   mem = mean(mem)) |> 
  mutate(method_clean = method_labs[method])

## L1 accuracy

l1_plot <- l1_cdf_summary |> 
  ggplot(aes(x = method, y = mean_l1)) +
  geom_boxplot() +
  lims(y = c(0.45, 1)) +
  scale_x_discrete(labels = method_labs) +
  labs(x = "", y = "Mean L1 accuracy") +
  theme_bw()

ggsave("plots/l1_plot.png", plot = l1_plot, dpi = 600, width = 24, height = 14, units = "cm")

## Time usage

time_plot <- other_cdf |> 
  ggplot(aes(x = method, y = time)) +
  geom_boxplot() +
  scale_x_discrete(labels = method_labs) +
  scale_y_log10() +
  labs(x = "", y = "Run time in seconds (log scale)") +
  theme_bw()

ggsave("plots/time_plot.png", plot = time_plot, dpi = 600, width = 24, height = 14, units = "cm")

## Memory usage

mem_plot <- other_cdf |> 
  ggplot(aes(x = method, y = mem)) +
  geom_boxplot() +
  scale_x_discrete(labels = method_labs) +
  scale_y_log10() +
  labs(x = "", y = "Memory usage in bytes (log scale)") +
  theme_bw()

ggsave("plots/mem_plot.png", plot = mem_plot, dpi = 600, width = 24, height = 14, units = "cm")

## Combined plot

combined_plot <- mean_cdf |> 
  ggplot(aes(x = time, y = mean_l1, size = log10(mem), label = method_clean)) +
  geom_point() +
  geom_text(size = 4, hjust = 0.5, vjust = -1.2) +
  lims(y = c(0.7, 1)) +
  scale_x_log10() +
  labs(x = "Run time in seconds (log scale)", 
       y = "Mean L1 accuracy",
       size = "Log base 10 of memory usage in bytes") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("plots/combined_plot.png", plot = combined_plot, dpi = 600, width = 24, height = 14, units = "cm")
