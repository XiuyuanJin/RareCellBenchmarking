## =========================
## Figure_4
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section4/figure")

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(patchwork)
library(purrr)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 4a ----
## =========================
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

prop_breaks <- c(0.00125,0.0025,0.005,0.01,0.02,0.03,0.04,0.05)
prop_labels <- c("0.125%","0.25%","0.5%","1%","2%","3%","4%","5%")

scenario_method_F1 <- data %>%
  group_by(Scenario,SampleID,Celltype,Cellname,PrecRare_num,PrecRare,Method) %>%
  summarise(
    median_F1 =
      if (all(is.na(F1))) {
        NA_real_
      } else {
        median(F1,na.rm = TRUE)
      },
    .groups = "drop"
  )

scenario_silhouette <- data %>%
  group_by(Scenario,SampleID,Celltype,Cellname,PrecRare_num,PrecRare) %>%
  summarise(
    Rare_cell_silhouette =
      if (all(is.na(Rare_cell_silhouette))) {
        NA_real_
      } else {
        median(Rare_cell_silhouette,na.rm = TRUE)
      },
    .groups = "drop"
  )

scenario_data <- scenario_method_F1 %>%
  left_join(scenario_silhouette %>% dplyr::select(Scenario,Rare_cell_silhouette), by = "Scenario")

scenario_data$Method <- factor(scenario_data$Method, levels = method_order)

scenario_silhouette_group <- scenario_silhouette %>%
  filter(!is.na(Rare_cell_silhouette)) %>%
  group_by(PrecRare_num) %>%
  mutate(silhouette_tertile = ntile(Rare_cell_silhouette,3),
         Silhouette_group = case_when(
           silhouette_tertile == 1 ~ "Low",
           silhouette_tertile == 2 ~ "Intermediate",
           silhouette_tertile == 3 ~ "High")
  ) %>%
  ungroup() %>%
  mutate(Silhouette_group = factor(Silhouette_group,levels = c("Low","Intermediate","High")))

scenario_data2 <- scenario_method_F1 %>%
  left_join(scenario_silhouette_group %>% dplyr::select(Scenario, Rare_cell_silhouette, Silhouette_group), by = "Scenario")

scenario_data2$Method <- factor(scenario_data2$Method, levels = method_order)

plot_a_data <- scenario_data2 %>%
  filter(
    !is.na(median_F1),
    !is.na(Silhouette_group)
  ) %>%
  group_by(Method, PrecRare_num, PrecRare, Silhouette_group) %>%
  summarise(
    F1_median =median(median_F1, na.rm = TRUE),
    F1_Q1 = quantile(median_F1, 0.25, na.rm = TRUE),
    F1_Q3 = quantile(median_F1, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


Silhouette_color <- c("Low" = "#C96E2D","Intermediate" = "#EFCB68","High" = "#1C4E73")

p1 <- ggplot(plot_a_data, aes(x = PrecRare,y = F1_median,color = Silhouette_group,fill = Silhouette_group,group = Silhouette_group)) +
  geom_ribbon(aes(ymin = F1_Q1,ymax = F1_Q3),alpha = 0.10,color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ Method,ncol = 4) +
  scale_color_manual(values = Silhouette_color) +
  scale_fill_manual(values = Silhouette_color) +
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25)) +
  labs(x = "Rare-cell proportion", y = "Median F1 score") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 11, color = "black"),
        legend.position = "none",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid = element_blank())

p1
ggsave(p1, filename = "Fig4.silhouette.lineplot.pdf",width= 12, height= 5.5, units='in')


## =========================
## Figure S6 ----
## =========================
scenario_method_median <- data %>%
  group_by(Scenario,SampleID,Celltype,Cellname,PrecRare_num,PrecRare,Method) %>%
  summarise(
    median_F1 = if (
      all(is.na(F1))
    ) {
      NA_real_
    } else {
      median(F1, na.rm = TRUE)
    },
    
    Rare_cell_silhouette = if (
      all(is.na(Rare_cell_silhouette))
    ) {
      NA_real_
    } else {
      median(
        Rare_cell_silhouette,
        na.rm = TRUE
      )
    },
    
    .groups = "drop"
  ) %>%
  mutate(Method = factor(Method,levels = method_order),
         PrecRare = factor(PrecRare,levels = prop_labels))


cor_stats <- scenario_method_median %>%
  filter(
    !is.na(median_F1),
    !is.na(Rare_cell_silhouette)
  ) %>%
  group_by(
    Method,
    PrecRare
  ) %>%
  group_modify(
    ~ {
      
      if (
        nrow(.x) < 3 ||
        n_distinct(.x$Rare_cell_silhouette) < 2 ||
        n_distinct(.x$median_F1) < 2
      ) {
        
        tibble(
          rho = NA_real_,
          p_value = NA_real_,
          n = nrow(.x)
        )
        
      } else {
        
        tmp <- suppressWarnings(
          cor.test(
            .x$Rare_cell_silhouette,
            .x$median_F1,
            method = "spearman",
            exact = FALSE
          )
        )
        
        tibble(
          rho = unname(tmp$estimate),
          p_value = tmp$p.value,
          n = nrow(.x)
        )
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    label = ifelse(
      is.na(rho),
      "",
      paste0(
        "R = ",
        sprintf("%.2f", rho),
        "\np = ",
        formatC(
          p_value,
          format = "g",
          digits = 3
        )
      )
    )
  )

plot_df <- scenario_method_median %>%
  filter(
    !is.na(median_F1),
    !is.na(Rare_cell_silhouette)
  )

p2 <- ggplot(plot_df,aes(x = Rare_cell_silhouette,y = median_F1,color = Method)) +
  geom_point(size = 0.75,alpha = 0.65) +
  geom_smooth(method = "lm",se = TRUE, color = "#2C3E50", fill = "#BFC9D4",linewidth = 0.8,alpha = 0.25) +
  geom_text(data = cor_stats,aes(x = -Inf,y = Inf,label = label),
            inherit.aes = FALSE,
            hjust = -0.10,
            vjust = 1.25,
            size = 2.7,
            fontface = "italic") +
  facet_grid(Method ~ PrecRare) +
  scale_color_manual(values = Methodcolor) +
  scale_x_continuous(
    limits = c(-0.45, 0.80),
    breaks = c(-0.25,0,0.25,0.5,0.75),
    expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0,0.5,1),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Rare-cell silhouette", y = "Median F1 score") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white",color = "black",linewidth = 0.4),
        strip.text.x = element_text(size = 9,face = "bold"),
        strip.text.y = element_text(size = 9,face = "bold"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 7),
        panel.spacing = unit(0.08,"cm"))
p2
ggsave(p2, filename = "FigS6.silhouette.all.pdf",width= 13, height= 9, units='in')


## =========================
## Figure 4b ----
## =========================
scenario_method_precision <- data %>%
  group_by(Scenario, SampleID, Celltype, Cellname, PrecRare_num, PrecRare, Method) %>%
  summarise(
    median_precision = if (all(is.na(precision))) NA_real_ else median(precision, na.rm = TRUE),
    .groups = "drop"
  )

precision_silhouette <- scenario_method_precision %>%
  left_join(
    scenario_silhouette_group %>%
      dplyr::select(Scenario, PrecRare_num, PrecRare, Rare_cell_silhouette, Silhouette_group),
    by = c("Scenario", "PrecRare_num", "PrecRare")
  )

precision_summary <- precision_silhouette %>%
  filter(!is.na(median_precision), !is.na(Silhouette_group)) %>%
  group_by(Method, PrecRare_num, PrecRare, Silhouette_group) %>%
  summarise(median_precision = median(median_precision, na.rm = TRUE), .groups = "drop")

precision_summary <- precision_summary %>%
  mutate(Method = factor(Method, levels = method_order),
         Silhouette_group = factor(Silhouette_group, levels = c("Low", "Intermediate", "High")))

p3 <- ggplot(precision_summary,aes(x = Method, y = median_precision, fill = Silhouette_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.75,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_point(
    aes(group = Silhouette_group),
    position = position_jitterdodge(
      jitter.width = 0.05,
      dodge.width = 0.75
    ),
    size = 1.4,
    alpha = 0.65,
    shape = 21
  ) +
  scale_fill_manual(values = Silhouette_color) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(x = "", y = "Median precision", fill = "Rare-cell silhouette") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9,color = "black"),
    axis.text.y = element_text(size = 9,color = "black"),
    axis.title.y = element_text(size = 11,color = "black"),
    legend.position = "none",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

p3
ggsave(p3, filename = "Fig4b.precision.box.pdf",width= 11, height = 3, units='in')



## Precision ----
scenario_method_precision <- data %>%
  group_by(Scenario, SampleID, Celltype, Cellname, PrecRare_num, PrecRare, Method) %>%
  summarise(
    median_precision = if (all(is.na(precision))) NA_real_ else median(precision, na.rm = TRUE),
    .groups = "drop"
  )

precision_data <- scenario_method_precision %>%
  left_join(
    scenario_silhouette_group %>% dplyr::select(Scenario, Rare_cell_silhouette, Silhouette_group),
    by = "Scenario"
  ) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    PrecRare = factor(PrecRare, levels = prop_labels)
  )

plot_precision_data <- precision_data %>%
  filter(!is.na(median_precision), !is.na(Silhouette_group)) %>%
  group_by(Method, PrecRare_num, PrecRare, Silhouette_group) %>%
  summarise(
    precision_median = median(median_precision, na.rm = TRUE),
    precision_Q1 = quantile(median_precision, 0.25, na.rm = TRUE),
    precision_Q3 = quantile(median_precision, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_precision <- ggplot(
  plot_precision_data,
  aes(x = PrecRare, y = precision_median, color = Silhouette_group,
      fill = Silhouette_group, group = Silhouette_group)
) +
  geom_ribbon(aes(ymin = precision_Q1, ymax = precision_Q3), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~Method, ncol = 6) +
  scale_color_manual(values = Silhouette_color) +
  scale_fill_manual(values = Silhouette_color) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Rare-cell proportion", y = "Median precision") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, color = "black"),
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid = element_blank()
  )

p_precision

ggsave("FigS7a.silhouette.precision.lineplot.pdf",p_precision, width = 13, height = 5, units = "in")



## Recall ----
scenario_method_recall <- data %>%
  group_by(Scenario, SampleID, Celltype, Cellname, PrecRare_num, PrecRare, Method) %>%
  summarise(
    median_recall = if (all(is.na(recall))) NA_real_ else median(recall, na.rm = TRUE),
    .groups = "drop"
  )

recall_data <- scenario_method_recall %>%
  left_join(
    scenario_silhouette_group %>% dplyr::select(Scenario, Rare_cell_silhouette, Silhouette_group),
    by = "Scenario"
  ) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    PrecRare = factor(PrecRare, levels = prop_labels)
  )

plot_recall_data <- recall_data %>%
  filter(!is.na(median_recall), !is.na(Silhouette_group)) %>%
  group_by(Method, PrecRare_num, PrecRare, Silhouette_group) %>%
  summarise(
    recall_median = median(median_recall, na.rm = TRUE),
    recall_Q1 = quantile(median_recall, 0.25, na.rm = TRUE),
    recall_Q3 = quantile(median_recall, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p_recall <- ggplot(
  plot_recall_data,
  aes(x = PrecRare, y = recall_median, color = Silhouette_group,
      fill = Silhouette_group, group = Silhouette_group)
) +
  geom_ribbon(aes(ymin = recall_Q1, ymax = recall_Q3), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~Method, ncol = 6) +
  scale_color_manual(values = Silhouette_color) +
  scale_fill_manual(values = Silhouette_color) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Rare-cell proportion", y = "Median recall") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, color = "black"),
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid = element_blank()
  )

p_recall

ggsave("FigS7b.silhouette.recall.lineplot.pdf", p_recall, width = 13, height = 5, units = "in")
