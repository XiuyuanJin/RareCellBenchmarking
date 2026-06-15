## =========================
## Figure_5
## =========================

library(dplyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/HLCA.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 5a ----
## =========================
df_line <- anno_data %>%
  group_by(method, ann_level) %>%
  dplyr::summarise(F1 = median(F1, na.rm = TRUE),
                   MCC = median(MCC, na.rm = TRUE),
                   recall = median(recall, na.rm = TRUE),
                   precision = median(precision, na.rm = TRUE))

p1 <- ggplot(df_line, aes(x = ann_level, y = F1, group = method, color = method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = Methodcolor) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  labs(x = "", y = "Median F1 score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p1


## =========================
## Figure 5b ----
## =========================
p2 <- ggplot(df_line, aes(x = ann_level, y = precision, group = method, color = method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = Methodcolor) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  labs(x = "", y = "Median precision") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p2


## =========================
## Figure 5c ----
## =========================
p3 <- ggplot(df_line, aes(x = ann_level, y = recall, group = method, color = method)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = Methodcolor) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic() +
  labs(x = "", y = "Median recall") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p3


## =========================
## Figure S7_1 ----
## =========================
df_summary <- anno_data %>%
  filter(ann_level == "ann_level_2") %>%
  group_by(rare_cell_type, method) %>%
  summarise(median_F1 = median(F1, na.rm = TRUE))


p4 <- ggplot(df_summary, aes(x = method, y = median_F1, fill = method)) +
  geom_col(color = "black", linewidth = 0.1) +
  facet_wrap(~rare_cell_type, nrow = 1) +
  scale_fill_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = "", y = "F1-score", title = "ann_level_2") +
  theme(strip.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p4


## =========================
## Figure S7_2 ----
## =========================
df_summary <- anno_data %>%
  filter(ann_level == "ann_level_3") %>%
  group_by(rare_cell_type, method) %>%
  summarise(median_F1 = median(F1, na.rm = TRUE))

p5 <- ggplot(df_summary, aes(x = method, y = median_F1, fill = method)) +
  geom_col(color = "black", linewidth = 0.1) +
  facet_wrap(~rare_cell_type, nrow = 4) +
  scale_fill_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = "", y = "F1-score", title = "ann_level_3") +
  theme(strip.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p5


## =========================
## Figure S8 ----
## =========================
df_summary <- anno_data %>%
  filter(ann_level == "ann_level_4") %>%
  group_by(rare_cell_type, method) %>%
  summarise(median_F1 = median(F1, na.rm = TRUE))

p6 <- ggplot(df_summary, aes(x = method, y = median_F1, fill = method)) +
  geom_col(color = "black", linewidth = 0.1) +
  facet_wrap(~rare_cell_type, nrow = 6) +
  scale_fill_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = "", y = "F1-score", title = "ann_level_4") +
  theme(strip.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p6

