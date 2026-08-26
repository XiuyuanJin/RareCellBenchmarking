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
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(lmerTest)
library(emmeans)
library(tidyverse)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")

prop_breaks <- c(0.00125,0.0025,0.005,0.01,0.02,0.03,0.04,0.05)
prop_labels <- c("0.125%","0.25%","0.5%","1%","2%","3%","4%","5%")


## =========================
## Figure 4a ----
## =========================
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


cor_heatmap <- scenario_data %>%
  filter(
    !is.na(median_F1),
    !is.na(Rare_cell_silhouette)
  ) %>%
  group_by(Method, PrecRare_num) %>%
  summarise(
    n = n(),
    rho = if (
      n() >= 3 &&
      n_distinct(median_F1) > 1 &&
      n_distinct(Rare_cell_silhouette) > 1
    ) {
      cor(
        median_F1,
        Rare_cell_silhouette,
        method = "spearman"
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    PrecRare = factor(
      PrecRare_num,
      levels = prop_breaks,
      labels = prop_labels
    ),
    label = ifelse(
      is.na(rho),
      "",
      sprintf("%.2f", rho)
    )
  )


## rho matrix
rho_mat <- cor_heatmap %>%
  dplyr::select(Method, PrecRare_num, rho) %>%
  mutate(Method = as.character(Method)) %>%
  pivot_wider(names_from = PrecRare_num, values_from = rho) %>%
  column_to_rownames("Method") %>%
  as.matrix()

rho_mat <- rho_mat[, as.character(prop_breaks), drop = FALSE]
colnames(rho_mat) <- prop_labels

## label matrix
label_mat <- cor_heatmap %>%
  dplyr::select(Method, PrecRare_num, label) %>%
  mutate(Method = as.character(Method)) %>%
  pivot_wider(names_from = PrecRare_num, values_from = label) %>%
  column_to_rownames("Method") %>%
  as.matrix()

label_mat <- label_mat[, as.character(prop_breaks), drop = FALSE]
colnames(label_mat) <- prop_labels

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "white", "#B2182B")
)

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

ht <- Heatmap(
  rho_mat,
  name = "Spearman's \u03c1",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = method_order,
  
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  row_names_side = "left",
  column_names_side = "bottom",
  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  
  rect_gp = gpar(col = "white", lwd = 1),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      label_mat[i, j],
      x, y,
      gp = gpar(fontsize = 9, col = "black")
    )
  },
  
  heatmap_legend_param = list(
    title = "Spearman correlation coefficient (\u03c1)",
    at = c(-1, -0.5, 0, 0.5, 1),
    direction = "horizontal",
    legend_width = unit(5, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

pdf("Fig4a.Spearman_heatmap.pdf", width = 6, height = 5.5)
draw(ht, heatmap_legend_side = "top")
dev.off()


## =========================
## Figure 4b ----
## =========================
scissors_box <- scenario_data %>%
  filter(Method == "SCISSORS",
         PrecRare_num %in% c(0.00125, 0.05),
         !is.na(median_F1),
         !is.na(Rare_cell_silhouette)) %>%
  group_by(PrecRare_num) %>%
  mutate(Silhouette_group = ifelse(Rare_cell_silhouette <= median(Rare_cell_silhouette, na.rm = TRUE),
                                   "Low", "High")) %>%
  ungroup() %>%
  mutate(PrecRare = factor(PrecRare_num,levels = c(0.00125, 0.05),labels = c("0.125%", "5%")),
         Silhouette_group = factor(Silhouette_group, levels = c("Low", "High")))

scissors_summary <- scissors_box %>%
  group_by(PrecRare, Silhouette_group) %>%
  summarise(median_F1 = median(median_F1, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Silhouette_group, values_from = median_F1) %>%
  mutate(Delta_F1 = `High` - `Low`)

scissors_summary

scissors_label <- scissors_box %>%
  group_by(PrecRare) %>%
  summarise(y = max(median_F1, na.rm = TRUE) + 0.08, .groups = "drop") %>%
  left_join(scissors_summary %>% dplyr::select(PrecRare, Delta_F1), by = "PrecRare") %>%
  mutate(label = paste0("\u0394F1 = ", sprintf("%.2f", Delta_F1)))


p1 <- ggplot(scissors_box, aes(x = PrecRare, y = median_F1, fill = Silhouette_group)) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.55, alpha = 0.7) +
  geom_point(aes(color = Silhouette_group),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             size = 1.5, alpha = 0.8) +
  geom_text(data = scissors_label,
            aes(x = PrecRare, y = y, label = label),
            inherit.aes = FALSE,
            size = 4,
            color = "black") +
  scale_fill_manual(values = c("Low" = "#A9CFE8", "High" = "#18486B")) +
  scale_color_manual(values = c("Low" = "#A9CFE8", "High" = "#18486B")) +
  labs(x = "Rare-cell proportion", y = "F1 score", fill = NULL, color = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.position = "top")
p1
ggsave(p1, filename = "Fig4b.scissors.box.pdf",width= 4.5, height = 4.5, units='in')


## =========================
## Figure 4c ----
## =========================
methods_interaction <- setdiff(names(Methodcolor), "EDGE")

interaction_data <- scenario_data %>%
  filter(
    Method %in% methods_interaction,
    PrecRare_num %in% c(0.00125, 0.05),
    !is.na(median_F1),
    !is.na(Rare_cell_silhouette)
  ) %>%
  mutate(SampleID_Celltype = paste(SampleID, Celltype, Cellname, sep = "_"),
         Prop_group = factor(PrecRare_num, levels = c(0.00125, 0.05), labels = c("0.125%", "5%")))


interaction_result <- lapply(methods_interaction, function(method_i) {
  
  df_i <- interaction_data %>% filter(Method == method_i)
  
  fit <- lm(median_F1 ~ Rare_cell_silhouette * Prop_group, data = df_i)
  
  coef_table <- summary(fit)$coefficients
  
  interaction_term <- grep(
    "Rare_cell_silhouette:Prop_group",
    rownames(coef_table),
    value = TRUE
  )
  
  trends <- emtrends(
    fit,
    ~ Prop_group,
    var = "Rare_cell_silhouette"
  ) %>%
    as.data.frame()
  
  data.frame(
    Method = method_i,
    Slope_0.125 = trends$Rare_cell_silhouette.trend[trends$Prop_group == "0.125%"],
    Slope_5 = trends$Rare_cell_silhouette.trend[trends$Prop_group == "5%"],
    Gamma = coef_table[interaction_term, "Estimate"],
    P = coef_table[interaction_term, "Pr(>|t|)"]
  )
})


interaction_result <- bind_rows(interaction_result) %>%
  mutate(Significant = P < 0.05)

interaction_result

volcano_data <- interaction_result %>%
  mutate(
    Direction = case_when(
      P < 0.05 & Gamma > 0 ~ "Positive interaction",
      P < 0.05 & Gamma < 0 ~ "Negative interaction",
      TRUE ~ "Not significant"
    ),
    neglog10P = -log10(P)
  )

p2 <- ggplot(volcano_data, aes(x = Gamma, y = neglog10P)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.5) +
  geom_point(aes(color = Direction), size = 3.5) +
  ggrepel::geom_text_repel(
    aes(label = Method),
    size = 3.5,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Positive interaction" = "#B2182B",
      "Negative interaction" = "#2166AC",
      "Not significant" = "grey65"
    )
  ) +
  labs(x = expression("Interaction coefficient (" * gamma * ")"),
       y = expression(-log[10]("P"))) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 12, color = "black"))
    

p2
ggsave(p2, filename = "Fig4c.volcano.gamma.pdf",width= 4.5, height = 4.5, units='in')


## =========================
## Figure 4d_1
## =========================
scenario_method_median <- data %>%
  group_by(Scenario, SampleID, Celltype, Cellname, PrecRare_num, PrecRare, Method) %>%
  summarise(
    median_precision = if(all(is.na(precision))) NA_real_ else median(precision, na.rm = TRUE),
    Rare_cell_silhouette = if(all(is.na(Rare_cell_silhouette))) NA_real_ else median(Rare_cell_silhouette, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    PrecRare = factor(PrecRare, levels = prop_labels)
  )

cor_stats <- scenario_method_median %>%
  filter(!is.na(median_precision), !is.na(Rare_cell_silhouette)) %>%
  group_by(Method, PrecRare) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || n_distinct(.x$Rare_cell_silhouette) < 2 || n_distinct(.x$median_precision) < 2) {
      tibble(rho = NA_real_, p_value = NA_real_, n = nrow(.x))
    } else {
      tmp <- suppressWarnings(cor.test(.x$Rare_cell_silhouette, .x$median_precision,
                                       method = "spearman", exact = FALSE))
      tibble(rho = unname(tmp$estimate), p_value = tmp$p.value, n = nrow(.x))
    }
  }) %>%
  ungroup() %>%
  mutate(
    label = ifelse(
      is.na(rho), "",
      paste0("R = ", sprintf("%.2f", rho),
             "\np = ", formatC(p_value, format = "g", digits = 3))
    )
  )


plot_df_scCAD <- scenario_method_median %>%
  filter(Method == "scCAD",
         !is.na(median_precision),
         !is.na(Rare_cell_silhouette))

cor_stats_scCAD <- cor_stats %>% filter(Method == "scCAD")

p3 <- ggplot(plot_df_scCAD, aes(x = Rare_cell_silhouette, y = median_precision)) +
  geom_point(size = 1.1, alpha = 0.65, color = Methodcolor["scCAD"]) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2C3E50", fill = "#BFC9D4",
              linewidth = 0.8, alpha = 0.25) +
  geom_text(
    data = cor_stats_scCAD,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.10, vjust = 1.25,
    size = 3, fontface = "italic"
  ) +
  facet_grid(. ~ PrecRare) +
  scale_x_continuous(
    limits = c(-0.45, 0.80),
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(x = "Rare-cell silhouette", y = "Median precision") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    panel.spacing = unit(0.08, "cm"))

p3
ggsave(p3, filename = "Fig4d.scCAD.silhouette.precision.pdf",width= 13, height= 2.2, units='in')


## =========================
## Figure 4d_2
## =========================
scenario_method_median <- data %>%
  group_by(Scenario, SampleID, Celltype, Cellname, PrecRare_num, PrecRare, Method) %>%
  summarise(
    median_recall = if(all(is.na(recall))) NA_real_ else median(recall, na.rm = TRUE),
    Rare_cell_silhouette = if(all(is.na(Rare_cell_silhouette))) NA_real_ else median(Rare_cell_silhouette, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    PrecRare = factor(PrecRare, levels = prop_labels)
  )


cor_stats <- scenario_method_median %>%
  filter(!is.na(median_recall), !is.na(Rare_cell_silhouette)) %>%
  group_by(Method, PrecRare) %>%
  group_modify(~ {
    if(nrow(.x) < 3 || n_distinct(.x$Rare_cell_silhouette) < 2 || n_distinct(.x$median_recall) < 2) {
      tibble(rho = NA_real_, p_value = NA_real_, n = nrow(.x))
    } else {
      tmp <- suppressWarnings(cor.test(.x$Rare_cell_silhouette, .x$median_recall,
                                       method = "spearman", exact = FALSE))
      tibble(rho = unname(tmp$estimate), p_value = tmp$p.value, n = nrow(.x))
    }
  }) %>%
  ungroup() %>%
  mutate(
    label = ifelse(
      is.na(rho), "",
      paste0("R = ", sprintf("%.2f", rho),
             "\np = ", formatC(p_value, format = "g", digits = 3))
    )
  )


plot_df_scCAD <- scenario_method_median %>%
  filter(Method == "scCAD",
         !is.na(median_recall),
         !is.na(Rare_cell_silhouette))

cor_stats_scCAD <- cor_stats %>% filter(Method == "scCAD")

p4 <- ggplot(plot_df_scCAD, aes(x = Rare_cell_silhouette, y = median_recall)) +
  geom_point(size = 1.1, alpha = 0.65, color = Methodcolor["scCAD"]) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#2C3E50", fill = "#BFC9D4",
              linewidth = 0.8, alpha = 0.25) +
  geom_text(
    data = cor_stats_scCAD,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.10, vjust = 1.25,
    size = 3, fontface = "italic"
  ) +
  facet_grid(. ~ PrecRare) +
  scale_x_continuous(
    limits = c(-0.45, 0.80),
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(x = "Rare-cell silhouette", y = "Median recall") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
    strip.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    panel.spacing = unit(0.08, "cm"))

p4
ggsave(p4, filename = "Fig4d.scCAD.silhouette.recall.pdf",width= 13, height= 2.2, units='in')


