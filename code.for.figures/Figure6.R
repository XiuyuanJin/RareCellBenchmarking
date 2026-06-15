## =========================
## Figure_6
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")

## =========================
## Figure 6a ----
## =========================
df_summary_sample <- data %>%
  group_by(Method, PrecRare, SampleID) %>%
  dplyr::summarise(median_Time = median(Time, na.rm = TRUE))

p1 <- ggplot(df_summary_sample, aes(x = Method, y = median_Time, color = Method)) +
  geom_boxplot()+
  scale_color_manual(values = Methodcolor) +     
  theme_classic() +
  labs(x = "", y = "Running time (mins)", color = "Method") +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(color = "black",size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size = 10),
        axis.text.y = element_text(color = "black",size = 10))
p1   


## =========================
## Figure 6b ----
## =========================
df_summary <- data %>%
  mutate(Cell_Nbin = cut(Celltype_Ncell,
                         breaks = c(0, 2000, 3500, 5000, 7000, 9000, 11500, 15000),
                         labels = c("0-2k","2-3.5k","3.5-5k","5-7k","7-9k","9-11.5k","11.5-15k")))

df_summary_bin <- df_summary %>%
  group_by(Method, Cell_Nbin) %>%
  dplyr::summarise(Celltype_Ncell = median(Celltype_Ncell, na.rm = TRUE),
                   Time = median(Time, na.rm = TRUE))

p2 <- ggplot(df_summary_bin, aes(x = Cell_Nbin, y = Time, group = Method, color = Method)) +
  geom_line(aes(group = Method), size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = "Dataset size", y = "Running time (mins)", color = "Method") +
  theme(strip.background = element_blank(),
        legend.position = "right",
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
p2


## =========================
## Figure S9 ----
## =========================
df_summary <- data %>%
  group_by(Method, PrecRare, SampleID, Celltype_Ncell) %>%
  dplyr::summarise(median_Time = median(Time, na.rm = TRUE))

p3 <- ggplot(df_summary, aes(x = Celltype_Ncell, y = median_Time, color = Method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~SampleID, nrow = 3, scales = "free") +           
  scale_color_manual(values = Methodcolor) +
  labs(x = "Dataset size", y = "Running time (mins)", color = "Method")+          
  theme_classic() +
  theme(strip.background = element_blank(),       
        strip.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))
p3


## =========================
## Figure S10 ----
## =========================
df_summary_sample <- data %>%
  group_by(Method, PrecRare, SampleID) %>%
  dplyr::summarise(median_Time = median(Time, na.rm = TRUE))

p4 <- ggplot(df_summary_sample, aes(x = PrecRare, y = median_Time, color = Method, group = Method)) +
  geom_line(size = 0.8) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = Methodcolor) +
  facet_wrap(~SampleID, nrow = 3) +           
  theme_classic() +
  labs(x = "Rare cell proportion", y = "Running time (mins)", color = "Method") +
  theme(strip.background = element_blank(),       
        strip.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))
p4

