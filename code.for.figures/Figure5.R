## =========================
## Figure_5
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section5/figure")

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(ggrepel)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 5a
## =========================
time_summary <- data %>%
  group_by(Method, SampleID, Cellname, PrecRare) %>%
  summarise(Time = median(Time, na.rm = TRUE), .groups = "drop")

method_order_time <- time_summary %>%
  group_by(Method) %>%
  summarise(median_Time = median(Time, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_Time) %>%
  pull(Method)

time_summary$Method <- factor(time_summary$Method, levels = method_order_time)

adjacent_pairs <- tibble(
  group1 = head(method_order_time, -1),
  group2 = tail(method_order_time, -1)
)

wilcox_results <- map2_dfr(adjacent_pairs$group1, adjacent_pairs$group2, function(m1, m2) {
  paired_data <- inner_join(
    time_summary %>% filter(Method == m1) %>% dplyr::select(SampleID, Cellname, PrecRare, Time1 = Time),
    time_summary %>% filter(Method == m2) %>% dplyr::select(SampleID, Cellname, PrecRare, Time2 = Time),
    by = c("SampleID", "Cellname", "PrecRare")
  )
  
  test_result <- wilcox.test(paired_data$Time1, paired_data$Time2, paired = TRUE, exact = FALSE)
  tibble(group1 = m1, group2 = m2, n = nrow(paired_data), p = test_result$p.value)
}) %>%
  mutate(
    p.signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    y.position = max(time_summary$Time, na.rm = TRUE) * 1.35
  )

wilcox_results

wilcox_results <- wilcox_results %>%
  mutate(
    x1 = match(group1, method_order_time),
    x2 = match(group2, method_order_time),
    y.position = max(time_summary$Time, na.rm = TRUE) * 5,
    y.tip = y.position / 1.12,
    y.label = y.position * 1.08)


p1 <- ggplot(time_summary, aes(x = Method, y = Time, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.9, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8, linewidth = 0.4) +
  geom_segment(data = wilcox_results, aes(x = x1, xend = x2, y = y.position, yend = y.position),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = wilcox_results, aes(x = x1, xend = x1, y = y.position, yend = y.tip),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_segment(data = wilcox_results, aes(x = x2, xend = x2, y = y.position, yend = y.tip),
               inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = wilcox_results, aes(x = (x1 + x2) / 2, y = y.label, label = p.signif),
            inherit.aes = FALSE, size = 4) +
  scale_y_log10(
    breaks = 10^(-1:3),
    labels = scales::trans_format(
      "log10",
      scales::math_format(10^.x)
    ),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_manual(values = Methodcolor) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  labs(x = NULL, y = "Running time (minutes)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12))
p1
ggsave(p1, filename = "Fig5a.time.pdf",width= 8, height= 5, units='in')


## =========================
## Figure 5b
## =========================
time_cell_summary <- data %>% 
  group_by( Method, SampleID, Cellname, PrecRare, Total_cell_number ) %>%
  summarise( Time = mean(Time, na.rm = TRUE), .groups = "drop" )

time_cell_summary <- time_cell_summary %>%
  mutate(
    Cell_number_bin = cut(
      Total_cell_number,
      breaks = c(
        0,
        2000,
        4000,
        6000,
        9000,
        12000,
        Inf
      ),
      labels=c(
        "<2k",
        "2-4k",
        "4-6k",
        "6-9k",
        "9-12k",
        ">12k"
      ),
      right = FALSE
    )
  )

time_cell_summary$Cell_number_bin <- factor(
  time_cell_summary$Cell_number_bin,
  levels = c(
    "<2k",
    "2-4k",
    "4-6k",
    "6-9k",
    "9-12k",
    ">12k"
  )
)

time_bin_summary <- time_cell_summary %>%
  group_by(Method, Cell_number_bin) %>%
  summarise(Median_Time = median(Time), .groups = "drop")

p2 <- ggplot(time_bin_summary,aes(x = Cell_number_bin,y = Median_Time,color = Method,group = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = Methodcolor) +    
  theme_classic() +
  labs(x = "Number of cells", y = "Median runtime (min)") +
  theme(axis.title = element_text(color = "black",size = 14),
        axis.text = element_text(color = "black",size = 12))
p2
ggsave(p2, filename = "Fig5b.time.datasize.pdf",width= 6, height= 4, units='in')


## =========================
## Figure 5c
## =========================
time_scaling_data <- data %>%
  group_by(
    Method,
    SampleID,
    Cellname,
    PrecRare,
    Total_cell_number
  ) %>%
  summarise(
    Time = median(Time, na.rm = TRUE),
    .groups = "drop"
  )

runtime_scaling_sample <- time_scaling_data %>%
  group_by(Method) %>%
  filter(n_distinct(Total_cell_number) > 2) %>%
  do(
    tidy(
      lm(
        log10(Time)~
          log10(Total_cell_number),
        data=.
      )
    )
  ) %>%
  filter(term=="log10(Total_cell_number)")



beta_order <- runtime_scaling_sample %>%
  arrange(estimate) %>%
  pull(Method)

beta_df <- runtime_scaling_sample %>%
  mutate(Method = factor(Method, levels = rev(beta_order)),
         beta_label = sprintf("%.2f", estimate))

p3 <- ggplot(beta_df, aes(x = estimate, y = Method, color = Method)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = estimate, yend = Method), color = "grey70", linewidth = 0.6) +
  geom_point(size = 5) +
  geom_text(aes(label = beta_label), hjust = -0.5, size = 3.5, color = "black") +
  scale_color_manual(values = Methodcolor) +
  scale_x_continuous(limits = c(0, 2.6), breaks = c(0, 0.5, 1, 1.5, 2, 2.5),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression("Scaling coefficient (" * beta * ")"), y = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 11))
    
p3
ggsave(p3, filename = "Fig5c.time.datasize.pdf",width= 9, height= 4, units='in')

