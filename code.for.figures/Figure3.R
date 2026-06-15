## =========================
## Figure_3
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 3a ----
## =========================
data$SampleID_Celltype <- paste(data$SampleID, data$Celltype, sep = "_")

benchmark_prop <- data %>%
  group_by(SampleID_Celltype, SampleID, Celltype, Cellname,Method, PrecRare, PrecRare_num) %>%
  dplyr::summarise(
    F1        = median(F1, na.rm = TRUE),
    MCC       = median(MCC, na.rm = TRUE),
    precision = median(precision, na.rm = TRUE),
    recall    = median(recall, na.rm = TRUE),
    .groups = "drop")


p1 <- ggplot(benchmark_prop,aes(x = PrecRare, y = F1, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 4) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(x = "Rare-cell proportion",y = "F1 score") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 11, color = "black"))
p1


## =========================
## Figure 3b ----
## =========================
calc_slope <- function(df) {
  df <- df %>%
    dplyr::select(PrecRare_num, F1) %>%
    filter(!is.na(PrecRare_num), !is.na(F1)) %>%
    arrange(PrecRare_num)
  
  if (nrow(df) < 3 || n_distinct(df$PrecRare_num) < 3) {
    return(tibble(
      slope = NA_real_,
      r_squared = NA_real_,
      delta_low_high = NA_real_
    ))
  }
  
  x <- log10(df$PrecRare_num)
  y <- df$F1
  
  fit <- lm(y ~ x)
  
  tibble(
    slope = as.numeric(coef(fit)[2]),
    r_squared = summary(fit)$r.squared,
    delta_low_high = y[which.max(df$PrecRare_num)] - y[which.min(df$PrecRare_num)]
  )
}


sensitivity_F1 <- benchmark_prop %>%
  group_by(SampleID_Celltype, Method) %>%
  group_modify(~ calc_slope(.x)) %>%
  ungroup()

method_order_slope <- sensitivity_F1 %>%
  group_by(Method) %>%
  dplyr::summarise(median_slope = median(slope, na.rm = TRUE),
                   .groups = "drop") %>%
  arrange(desc(median_slope)) %>%
  pull(Method)

sensitivity_F1_plot <- sensitivity_F1 %>%
  mutate(Method = factor(Method, levels = method_order_slope))

p2 <- ggplot(sensitivity_F1_plot, aes(x = Method, y = slope, fill = Method)) +
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey45",linewidth = 0.5) +
  geom_boxplot(width = 0.7, alpha = 0.65, outlier.size = 0.5) +
  geom_point(aes(color = Method),alpha = 0.3,size = 0.5,position = position_jitter(width = 0.15)) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(x = "", y = "Slope of F1 vs log10(rare-cell proportion)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"))
p2


## =========================
## Figure S3 ----
## =========================
p3 <- ggplot(benchmark_prop,aes(x = PrecRare, y = precision, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 4) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(x = "Rare-cell proportion",y = "Precision") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 11, color = "black"))
p3


## =========================
## Figure S4 ----
## =========================
p4 <- ggplot(benchmark_prop,aes(x = PrecRare, y = recall, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 4) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(x = "Rare-cell proportion",y = "Recall") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 11, color = "black"))
p4
