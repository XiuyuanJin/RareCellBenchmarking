## =========================
## Figure_3
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section3/figure")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 3a-c ----
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
  "scCAD",
  "CellSIUS",
  "SCISSORS",
  "CIARA",
  "SCA",
  "GiniClust3",
  "RaceID3",
  "RaceID2",
  "EDGE"
)

benchmark_prop$Method <- factor(benchmark_prop$Method,levels = method_order)

p1 <- ggplot(benchmark_prop, aes(x = PrecRare, y = F1, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-F1 score") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p1
ggsave(p1, filename = "Fig3a.line.prop.F1.pdf",width= 6, height= 5, units='in')



p2 <- ggplot(benchmark_prop, aes(x = PrecRare, y = precision, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-precision") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p2
ggsave(p2, filename = "Fig3b.line.prop.precision.pdf",width= 6, height= 5, units='in')


p3 <- ggplot(benchmark_prop, aes(x = PrecRare, y = recall, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-recall") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p3
ggsave(p3, filename = "Fig3c.line.prop.recall.pdf",width= 6, height= 5, units='in')


## =========================
## Figure 3d-f ----
## =========================

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.scCAD.RData")

Methodcolor <- c(
  "scCAD(1%)" = "#08306B",
  "scCAD(2%)" = "#2171B5",
  "scCAD(3%)" = "#4292C6",
  "scCAD(4%)" = "#6BAED6",
  "scCAD(5%)" = "#C6DBEF"
)


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

p4 <- ggplot(benchmark_prop, aes(x = PrecRare, y = precision, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-precision") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p4
ggsave(p4, filename = "Fig3d.line.scCAD.precision.pdf",width= 6, height= 5, units='in')


p5 <- ggplot(benchmark_prop, aes(x = PrecRare, y = recall, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-recall") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p5
ggsave(p5, filename = "Fig3e.line.scCAD.recall.pdf",width= 6, height= 5, units='in')


p6 <- ggplot(benchmark_prop, aes(x = PrecRare, y = F1, color = Method, group = Method)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(values = Methodcolor) +
  ylim(c(0,1)) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "rs-F1 score") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p6
ggsave(p6, filename = "Fig3f.line.scCAD.F1.pdf",width= 6, height= 5, units='in')
