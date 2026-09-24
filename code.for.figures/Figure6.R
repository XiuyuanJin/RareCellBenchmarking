## =========================
## Figure_6
## =========================
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)

strategy_colors <- c("Single method" = "#0260A2","Two-method ensemble" = "#AD0F09","Three-method ensemble" = "#EEB542")

strategy_labels <- c(
  "aKNNO" = "aKNNO",
  "RareQ" = "RareQ",
  "CellSIUS" = "CellSIUS",
  
  "aKNNO_RareQ_union" =
    "aKNNO + RareQ\nunion",
  
  "aKNNO_RareQ_intersection" =
    "aKNNO + RareQ\nintersection",
  
  "aKNNO_CellSIUS_union" =
    "aKNNO + CellSIUS\nunion",
  
  "aKNNO_CellSIUS_intersection" =
    "aKNNO + CellSIUS\nintersection",
  
  "RareQ_CellSIUS_union" =
    "RareQ + CellSIUS\nunion",
  
  "RareQ_CellSIUS_intersection" =
    "RareQ + CellSIUS\nintersection",
  
  "three_union" =
    "Three-method\nunion",
  
  "three_majority" =
    "Three-method\nmajority",
  
  "three_intersection" =
    "Three-method\nintersection"
)



meta_data <- readRDS("/data/Figure6_metadata.rds")

samplename <- names(table(meta_data$sample_name))

base_dir <- "/data"
methods <- c("aKNNO", "RareQ", "CellSIUS")

read_one_sample <- function(sample_id) {
  
  finalclusters_aKNNO <- read.table(paste0(base_dir,"/Figure6_aKNNO/FetalLung_2024NC_Wong_",sample_id,"_aKNNO_finalcluster.txt"),
                                    sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  finalclusters_RareQ <- read.table(paste0(base_dir,"/Figure6_RareQ/FetalLung_2024NC_Wong_",sample_id,"_RareQ_finalcluster.txt"),
                                    sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  finalclusters_CellSIUS <- read.table(paste0(base_dir,"/Figure6_CellSIUS/FetalLung_2024NC_Wong_",sample_id,"_CellSIUS_finalcluster.txt"),
                                       sep = ",", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  
  cell_id <- rownames(finalclusters_aKNNO)
  
  df <- data.frame(sample = sample_id,
                   cell_id = cell_id,
                   aKNNO_cluster = as.character(finalclusters_aKNNO$cluster),
                   RareQ_cluster = as.character(finalclusters_RareQ$cluster),
                   CellSIUS_cluster = as.character(finalclusters_CellSIUS$cluster),
                   stringsAsFactors = FALSE)
  
  return(df)
}

cluster_wide <- map_dfr(samplename, read_one_sample)


## =========================
## Figure 6a ----
## =========================
rare_cut <- 0.01
min_cells <- 3

cluster_long <- cluster_wide %>%
  pivot_longer(
    cols = c(aKNNO_cluster, RareQ_cluster, CellSIUS_cluster),
    names_to = "method",
    values_to = "cluster"
  ) %>%
  mutate(
    method = sub("_cluster$", "", method),
    cluster = as.character(cluster)
  )

cluster_candidate <- cluster_long %>%
  group_by(sample, method) %>%
  mutate(total_cells = n_distinct(cell_id)) %>%
  group_by(sample, method, cluster) %>%
  mutate(
    cluster_cells = n_distinct(cell_id),
    cluster_prop = cluster_cells / dplyr::first(total_cells)
  ) %>%
  ungroup() %>%
  mutate(
    candidate_rare = cluster_cells >= min_cells &
      cluster_prop <= rare_cut
  )


ensemble_vote <- cluster_candidate %>%
  group_by(sample, cell_id) %>%
  summarise(aKNNO = any(method == "aKNNO" & candidate_rare, na.rm = TRUE),
            RareQ = any(method == "RareQ" & candidate_rare, na.rm = TRUE),
            CellSIUS = any(method == "CellSIUS" & candidate_rare, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(vote = as.integer(aKNNO) + as.integer(RareQ) + as.integer(CellSIUS),
         
         support_pattern = case_when(
           vote == 0 ~ "None",
           
           aKNNO & !RareQ & !CellSIUS ~ "aKNNO only",
           !aKNNO & RareQ & !CellSIUS ~ "RareQ only",
           !aKNNO & !RareQ & CellSIUS ~ "CellSIUS only",
           
           aKNNO & RareQ & !CellSIUS ~ "aKNNO + RareQ",
           aKNNO & !RareQ & CellSIUS ~ "aKNNO + CellSIUS",
           !aKNNO & RareQ & CellSIUS ~ "RareQ + CellSIUS",
           
           aKNNO & RareQ & CellSIUS ~ "All three"
         ),
         
         support_level = case_when(
           vote == 0 ~ "Not detected",
           vote == 1 ~ "Single-method detection",
           vote == 2 ~ "Two-method consensus",
           vote == 3 ~ "Three-method consensus"
         ),
         
         aKNNO_RareQ_union = aKNNO | RareQ,
         aKNNO_RareQ_intersection = aKNNO & RareQ,
         aKNNO_CellSIUS_union = aKNNO | CellSIUS,
         aKNNO_CellSIUS_intersection = aKNNO & CellSIUS,
         RareQ_CellSIUS_union = RareQ | CellSIUS,
         RareQ_CellSIUS_intersection = RareQ & CellSIUS,
         
         three_union = aKNNO | RareQ | CellSIUS,
         three_majority = vote >= 2,
         three_intersection = aKNNO & RareQ & CellSIUS
  )


strategy_columns <- c(
  "aKNNO",
  "RareQ",
  "CellSIUS",
  
  "aKNNO_RareQ_union",
  "aKNNO_RareQ_intersection",
  "aKNNO_CellSIUS_union",
  "aKNNO_CellSIUS_intersection",
  "RareQ_CellSIUS_union",
  "RareQ_CellSIUS_intersection",
  
  "three_union",
  "three_majority",
  "three_intersection"
)


ensemble_long <- ensemble_vote %>%
  dplyr::select(
    sample,
    cell_id,
    vote,
    support_pattern,
    support_level,
    all_of(strategy_columns)
  ) %>%
  pivot_longer(
    cols = all_of(strategy_columns),
    names_to = "strategy",
    values_to = "predicted_rare"
  ) %>%
  mutate(
    strategy_type = case_when(
      strategy %in% c("aKNNO", "RareQ", "CellSIUS") ~
        "Single method",
      
      grepl("^three_", strategy) ~
        "Three-method ensemble",
      
      TRUE ~
        "Two-method ensemble"
    )
  )


min_cells <- 3

truth <- meta_data %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  transmute(
    cell_id = as.character(cell_id),
    sample = as.character(sample_name),
    all_cell_type = as.character(all_cell_type)
  ) %>%
  filter(!is.na(all_cell_type)) %>%
  group_by(sample, all_cell_type) %>%
  mutate(n_celltype = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(
    total_cells = n_distinct(cell_id),
    celltype_prop = n_celltype / total_cells,
    annotation_rare = celltype_prop <= rare_cut,
    
    evaluation_eligible =
      !(annotation_rare & n_celltype < min_cells),
    
    true_rare_cell =
      annotation_rare & n_celltype >= min_cells
  ) %>%
  ungroup()


evaluation_data <- ensemble_long %>%
  left_join(
    truth %>%
      dplyr::select(
        sample,
        cell_id,
        all_cell_type,
        n_celltype,
        celltype_prop,
        evaluation_eligible,
        true_rare_cell
      ),
    by = c("sample", "cell_id")
  )

evaluation_data <- evaluation_data %>%
  filter(evaluation_eligible)



evaluation_metrics <- evaluation_data %>%
  group_by(sample, strategy_type, strategy) %>%
  summarise(
    TP = sum(predicted_rare & true_rare_cell),
    FP = sum(predicted_rare & !true_rare_cell),
    FN = sum(!predicted_rare & true_rare_cell),
    TN = sum(!predicted_rare & !true_rare_cell),
    
    n_predicted = sum(predicted_rare),
    n_true_rare = sum(true_rare_cell),
    
    precision = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
    recall = ifelse(TP + FN > 0, TP / (TP + FN), NA_real_),
    
    F1 = ifelse(2 * TP + FP + FN > 0, 2 * TP / (2 * TP + FP + FN), NA_real_),
    F0.5 = ifelse(1.25 * TP + 0.25 * FN + FP > 0, 1.25 * TP / (1.25 * TP + 0.25 * FN + FP), NA_real_),
    
    .groups = "drop"
  )

make_metric_boxplot <- function(metric_name, add_wilcox = FALSE) {
  
  plot_data <- evaluation_metrics %>%
    transmute(sample, strategy_type, strategy, value = .data[[metric_name]]) %>%
    filter(!is.na(value))
  
  strategy_order_df <- plot_data %>%
    group_by(strategy) %>%
    summarise(median_value = median(value), q25_value = quantile(value, 0.25), .groups = "drop") %>%
    arrange(desc(median_value), desc(q25_value))
  
  ensemble_order <- strategy_order_df %>%
    filter(!strategy %in% c("aKNNO", "RareQ", "CellSIUS")) %>%
    pull(strategy)
  
  single_order <- strategy_order_df %>%
    filter(strategy %in% c("aKNNO", "RareQ", "CellSIUS")) %>%
    pull(strategy)
  
  strategy_order <- c(ensemble_order, single_order)
  plot_data <- plot_data %>% mutate(strategy = factor(strategy, levels = strategy_order))
  
  if (add_wilcox) {
    
    wilcox_all <- purrr::map_dfr(ensemble_order, function(s) {
      paired_data <- plot_data %>%
        mutate(strategy = as.character(strategy)) %>%
        filter(strategy %in% c(s, "aKNNO")) %>%
        dplyr::select(sample, strategy, value) %>%
        tidyr::pivot_wider(names_from = strategy, values_from = value, values_fn = median)
      
      keep <- complete.cases(paired_data[, c(s, "aKNNO")])
      test <- wilcox.test(paired_data[[s]][keep], paired_data[["aKNNO"]][keep],
                          paired = TRUE, exact = FALSE, alternative = "two.sided")
      
      tibble(
        strategy = s,
        n = sum(keep),
        median_difference = median(paired_data[[s]][keep] - paired_data[["aKNNO"]][keep]),
        p = test$p.value
      )
    }) %>%
      mutate(
        p.adj.signif = case_when(
          p < 0.001 ~ "***",
          p < 0.01  ~ "**",
          p < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
    
    annotation_strategies <- c(
      "three_intersection",
      "aKNNO_CellSIUS_intersection",
      "RareQ_CellSIUS_intersection",
      "aKNNO_RareQ_intersection",
      "three_majority"
    )
    
    missing_strategies <- setdiff(annotation_strategies, ensemble_order)
    
    wilcox_annotation <- wilcox_all %>%
      filter(strategy %in% annotation_strategies) %>%
      mutate(
        x1 = match(strategy, strategy_order),
        x2 = match("aKNNO", strategy_order),
        span = abs(x2 - x1)
      ) %>%
      arrange(span) %>%
      mutate(
        y.position = seq(1.05, 1.25, length.out = n()),
        y.tip = y.position - 0.02,
        y.label = y.position + 0.015
      )
    
    print(wilcox_all)
  }
  
  y_upper <- if (add_wilcox) 1.28 else 1
  
  p <- ggplot(plot_data, aes(x = strategy, y = value, fill = strategy_type)) +
    geom_point(aes(color = strategy_type), size = 0.8,
               position = position_jitter(width = 0.15)) +
    geom_boxplot(width = 0.7, alpha = 0.6) +
    scale_y_continuous(limits = c(0, y_upper), breaks = seq(0, 1, 0.25)) +
    scale_x_discrete(labels = strategy_labels) +
    scale_fill_manual(values = strategy_colors) +
    scale_color_manual(values = strategy_colors) +
    labs(x = NULL, y = metric_name) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, color = "black")
    )
  
  if (add_wilcox) {
    p <- p +
      geom_segment(data = wilcox_annotation,
                   aes(x = x1, xend = x2, y = y.position, yend = y.position),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = wilcox_annotation,
                   aes(x = x1, xend = x1, y = y.position, yend = y.tip),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = wilcox_annotation,
                   aes(x = x2, xend = x2, y = y.position, yend = y.tip),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(data = wilcox_annotation,
                aes(x = (x1 + x2) / 2, y = y.label, label = p.adj.signif),
                inherit.aes = FALSE, size = 4)
  }
  
  p
}


p1 <- make_metric_boxplot("precision", add_wilcox = TRUE)
p1
ggsave(p1, filename = "/results/Fig6a.precision.0.1.sig.pdf",width= 7, height= 4.5, units='in')



