## =========================
## Figure_3
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section3/figure")

library(dplyr)
library(tidyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")


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


method_order <- c(
  "aKNNO",
  "RareQ",
  "SCISSORS",
  "CIARA",
  
  "CellSIUS",
  "GiniClust3",
  "SCA",
  "scCAD",
  
  "RaceID2",
  "RaceID3",
  "EDGE"
)

benchmark_prop$Method <- factor(benchmark_prop$Method,levels = method_order)

p1 <- ggplot(benchmark_prop,aes(x = PrecRare, y = F1, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  geom_hline(yintercept = 0.75,linetype = 2,colour = "grey50") +
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
ggsave(p1, filename = "Fig3.prop.F1.pdf",width= 12, height= 5.5, units='in')


p2 <- ggplot(benchmark_prop,aes(x = PrecRare, y = precision, fill = Method)) +
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
p2
ggsave(p2, filename = "FigS4a.prop.precision.pdf",width= 12, height= 5.5, units='in')


p3 <- ggplot(benchmark_prop,aes(x = PrecRare, y = recall, fill = Method)) +
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
p3
ggsave(p3, filename = "FigS5a.prop.recall.pdf",width= 12, height= 5.5, units='in')


p4 <- ggplot(benchmark_prop, aes(x = PrecRare, y = precision, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "Median precision") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p4
ggsave(p4, filename = "Fig3c.line.prop.precision.pdf",width= 6, height= 5, units='in')


p5 <- ggplot(benchmark_prop, aes(x = PrecRare, y = recall, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "Median recall") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p5
ggsave(p5, filename = "Fig3d.line.prop.recall.pdf",width= 6, height= 5, units='in')


## =========================
## Figure S3b ----
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

p6 <- ggplot(sensitivity_F1_plot, aes(x = Method, y = slope, fill = Method)) +
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey45",linewidth = 0.5) +
  geom_boxplot(width = 0.7, alpha = 0.65, outlier.size = 0.5) +
  geom_point(aes(color = Method),alpha = 0.3,size = 0.5,position = position_jitter(width = 0.15)) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(x = "", y = "Slope of F1 vs log10\n(rare-cell proportion)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"))
p6
ggsave(p6, filename = "Fig3b.slope.F1.pdf",width= 7, height= 5, units='in')


## =========================
## Figure 3b ----
## =========================
load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.scCAD.RData")

Methodcolor <- c("scCAD(1%)" = "#08306B",
                 "scCAD(2%)" = "#2171B5",
                 "scCAD(3%)" = "#4292C6",
                 "scCAD(4%)" = "#6BAED6",
                 "scCAD(5%)" = "#C6DBEF")
  
data$SampleID_Celltype <- paste(data$SampleID, data$Celltype, sep = "_")

benchmark_prop <- data %>%
  group_by(SampleID_Celltype, SampleID, Celltype, Cellname,Method, PrecRare, PrecRare_num) %>%
  dplyr::summarise(
    F1        = median(F1, na.rm = TRUE),
    MCC       = median(MCC, na.rm = TRUE),
    precision = median(precision, na.rm = TRUE),
    recall    = median(recall, na.rm = TRUE),
    .groups = "drop")

method_order <- c("scCAD(1%)","scCAD(2%)","scCAD(3%)","scCAD(4%)","scCAD(5%)")
benchmark_prop$Method <- factor(benchmark_prop$Method,levels = method_order)

p1 <- ggplot(benchmark_prop, aes(x = PrecRare, y = F1, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "Median F1") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p1
ggsave(p1, filename = "Fig3b.line.scCAD.F1.pdf",width= 6, height= 5, units='in')


p2 <- ggplot(benchmark_prop,aes(x = PrecRare, y = F1, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 5) +
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
p2
ggsave(p2, filename = "FigS3a.scCAD.F1.pdf",width= 18, height= 3.5, units='in')


p3 <- ggplot(benchmark_prop,aes(x = PrecRare, y = precision, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 5) +
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
ggsave(p3, filename = "FigS4b.scCAD.precision.pdf",width= 18, height= 3.5, units='in')


p4 <- ggplot(benchmark_prop,aes(x = PrecRare, y = recall, fill = Method)) +
  geom_boxplot(width = 0.65,outlier.size = 0.3,linewidth = 0.25,alpha = 0.5) +
  stat_summary(aes(group = Method, color = Method),fun = median,geom = "line",linewidth = 0.8) +
  stat_summary(aes(color = Method),fun = median,geom = "point",size = 1) +
  facet_wrap(~ Method, ncol = 5) +
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
ggsave(p4, filename = "FigS5b.scCAD.recall.pdf",width= 18, height= 3.5, units='in')

