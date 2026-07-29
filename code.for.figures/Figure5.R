## =========================
## Figure_5
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section5/figure")

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")

## =========================
## Figure 5a ----
## =========================
time_summary <- data %>%
  group_by(
    Method,
    SampleID,
    Cellname,
    PrecRare
  ) %>%
  summarise(
    Time = median(Time),
    .groups = "drop"
  )

p1 <- ggplot(time_summary, aes(x = reorder(Method, Time, FUN = median), y = Time, fill = Method)) +
  geom_boxplot(width = 0.7, outlier.size = 0.8, alpha = 0.6) +
  geom_point(aes(fill = Method), shape = 21,  size = 0.7, stroke = 0.1) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
  scale_fill_manual(values = Methodcolor) +    
  geom_hline(yintercept = 60, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  theme_classic() +
  labs(x = NULL, y = "Running time (min)") +
  scale_color_manual(values = Methodcolor) +
  theme(legend.position = "none",
        axis.title.y = element_text(color = "black",size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))
p1
ggsave(p1, filename = "Fig5a.time.pdf",width= 7, height= 5, units='in')


## =========================
## Figure 5b ----
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
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  scale_color_manual(values = Methodcolor) +    
  theme_classic() +
  labs(x = "Number of cells", y = "Median runtime (min)") +
  theme(axis.title = element_text(color = "black",size = 14),
        axis.text = element_text(color = "black",size = 12))
p2
ggsave(p2, filename = "FigS5b.time.datasize.pdf",width= 6, height= 4, units='in')


## =========================
## Figure S8a ----
## =========================
time_prop_summary <- data %>% 
  group_by(
    Method,
    SampleID,
    Cellname,
    PrecRare
  ) %>%
  summarise(
    Time = mean(Time, na.rm = TRUE),
    .groups = "drop"
  )

time_prop_summary$PrecRare <- factor(time_prop_summary$PrecRare,levels = c("0.125%","0.25%","0.5%","1%","2%","3%","4%","5%"))

time_prop_bin_summary <- time_prop_summary %>%
  group_by(
    Method,
    PrecRare
  ) %>%
  summarise(
    Median_Time = median(Time, na.rm = TRUE),
    .groups = "drop"
  )


p3 <- ggplot(time_prop_bin_summary, aes(x = PrecRare, y = Median_Time, color = Method, group = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  scale_color_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = "Rare-cell proportion", y = "Median runtime (mins)") +
  theme(axis.title = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))
p3
ggsave(p3, filename = "FigS8a.time.prop.pdf",width= 6, height= 4, units='in')


## =========================
## Figure S8b ----
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
  group_by(
    Method,
    SampleID
  ) %>%
  filter(
    n_distinct(Total_cell_number)>2
  ) %>%
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
  group_by(Method) %>%
  summarise(
    median_beta = median(estimate,na.rm=TRUE)
  ) %>%
  arrange(median_beta) %>%
  pull(Method)


runtime_scaling_sample$Method <- factor(
  runtime_scaling_sample$Method,
  levels=beta_order
)

p4 <- ggplot(runtime_scaling_sample,aes(x = Method,y = estimate,fill = Method)) +
  geom_boxplot(width = 0.7, outlier.size = 0.8, alpha = 0.6) +
  geom_point(aes(fill = Method), shape = 21,  size = 0.7, stroke = 0.1) +
  scale_fill_manual(values = Methodcolor) +
  theme_classic() +
  labs(x = NULL, y = "Runtime scaling coefficient (β)") +
  geom_hline(yintercept = 1,linetype = "dashed",color = "grey30",linewidth = 0.3) +
  theme(legend.position = "none",
        axis.title.y = element_text(color = "black",size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1,color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))
p4
ggsave(p4, filename = "FigS8b.beta.pdf",width= 7, height= 5, units='in')

