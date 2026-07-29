## =========================
## Figure_6
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section6/figure")

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")


## =========================
## Figure 6a ----
## =========================
df_overall <- data %>%
  group_by(Method) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    index = "Overall performance",
    rank = rank(-value, ties.method = "min")
  )


df_F1 <- data %>%
  mutate(
    index = cut(
      PrecRare_num,
      breaks = c(0, 0.0025, 0.02, 0.05),
      labels = c(
        "Extremely rare\n(0.125–0.25%)",
        "Moderately rare\n(0.5–2%)",
        "Slightly rare\n(3–5%)"
      ),
      include.lowest = TRUE
    )
  ) %>%
  group_by(Method, index) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(index) %>%
  mutate(
    rank = rank(-value, ties.method = "min")
  ) %>%
  ungroup()


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

df_silhouette <- data %>%
  left_join(
    scenario_silhouette_group %>% dplyr::select(Scenario, Silhouette_group),
    by = "Scenario"
  ) %>%
  filter(!is.na(Silhouette_group)) %>%
  group_by(Method, Silhouette_group) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Silhouette_group) %>%
  mutate(
    rank = rank(-value, ties.method = "min"),
    index = case_when(
      Silhouette_group == "Low" ~ "Low silhouette",
      Silhouette_group == "Intermediate" ~ "Intermediate silhouette",
      Silhouette_group == "High" ~ "High silhouette"
    )
  ) %>%
  ungroup()



df_time <- data %>%
  group_by(Method) %>%
  dplyr::summarise(
    value = median(Time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    index = "Runtime efficiency",
    rank = rank(value, ties.method = "min")
  )


df_all <- bind_rows(
  df_overall,
  df_F1,
  df_silhouette,
  df_time
)

context_order <- c(
  "Overall performance",
  "Extremely rare\n(0.125–0.25%)",
  "Moderately rare\n(0.5–2%)",
  "Slightly rare\n(3–5%)",
  "Low silhouette",
  "Intermediate silhouette",
  "High silhouette",
  "Runtime efficiency"
)

df_all <- df_all %>%
  mutate(index = factor(index, levels = context_order))


method_order <- df_all %>%
  group_by(Method) %>%
  dplyr::summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_rank) %>%
  pull(Method)

df_all <- df_all %>%
  mutate(Method = factor(Method, levels = rev(method_order)))


p1 <- ggplot(df_all, aes(x = index, y = Method, fill = rank)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = rank), size = 3.2, color = "black") +
  scale_fill_gradientn(
    colours = c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7"),
    values = scales::rescale(c(1, 3, 5, 8, 11)),
    limits = c(1, 11),
    breaks = c(1, 3, 5, 7, 9, 11),
    name = "Rank") +
  labs(x = "", y = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 11),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 9, color = "black"))
p1
ggsave(p1, filename = "Fig6a.heatmap.pdf",width= 7, height= 4, units='in')


