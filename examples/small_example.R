library(phiper)
library(ggplot2)
library(magrittr)
library(dplyr)

# ----------- loading the data -------------------------------------------------
ps <- phip_convert(
  data_long_path = "../local-data/babies_phip_data.parquet",
  backend = "duckdb",
  peptide_library = TRUE,
  subject_id = "subject_id",
  peptide_id = "peptide_id",
  sample_id  = "sample_id",
  exist      = "exist",
  timepoint  = "timepoint_factor",
  fold_change= "fold_change",
  materialise_table = TRUE,
  auto_expand = FALSE,
  n_cores = 5
)

# get the current peptide library
get_peptide_library(ps)

# or alternatively with get_peptide_meta() --> it will download the library
# directly from my github

# ---------- plotting the enrichment counts ------------------------------------
plot_enrichment_counts(
  ps,
  group_cols = c("big_group", "timepoint"),
  group_interaction = TRUE,
  interaction_only = TRUE,
  annotation_size = 3
)

# ---------- alpha_diversity ---------------------------------------------------
alpha <- compute_alpha_diversity(ps,
                                 group_cols = c("big_group", "timepoint"),
                                 carry_cols = "subject_id",
                                 group_interaction = TRUE,
                                 interaction_only = TRUE)

# plot alpha diversity
plot_alpha_diversity(alpha,
                     metric = "richness",
                     group_col = "phip_interaction") +
  ylim(0, 2500) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# ---------- stability between kids at B and kids at M3 ------------------------
stab_recoded_tbl <- compute_repertoire_similarity(
  ps %>% filter(big_group == "kid_serum" & timepoint %in% c("T2", "T6")),
  group_col = "big_group",
  time_col = "timepoint",
  similarity = "kulczynski",
  mode = "all",
  dyade_col = "dyade_recoded",
  time_mode = "pairwise",
  time_pairing = c("same", "cross"),
  subject_mode = "all",
  drop_self_time = FALSE,
  overwrite = TRUE,
  verbose = TRUE
)

# plot the stabilities
plot_similarity_heatmap(stab_recoded_tbl,
                        groups = c("kid_serum", "kid_serum"),
                        times = c("T2", "T6"),
                        pairing = "cross",
                        dyads = "order") +
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1, # use -1 to invert
    limits = c(0, 1),
    oob = scales::squish,
    name = "Kulczynski similarity"
  )

# you can also plot them as boxplots; dyade = the same kid, you can group them
# by subject_id as well, its the same
sim <- stab_recoded_tbl %>%
  transmute(similarity,
            same = coalesce(dyad_left, "") == coalesce(dyad_right, "")) %>%
  collect()

ggplot(sim, aes(x = same, y = similarity)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0) +
  labs(x = NULL, y = "Kulczynski similarity") +
  theme_phip()

# ---------- compare prevalence using POP --------------------------------------
# lets say we want to compare different species for kids between B and M3
# minimal example because lightweight
POP_prev <- ph_prevalence_compare(
  ps %>% filter(big_group == "kid_serum" & timepoint %in% c("T2", "T6")),
  group_cols = "timepoint",
  rank_cols = "species",
  weight_mode = "peptide_count",
  paired = "subject_id"
)

# it is experimental and doesnt have any docs
scatter_static(POP_prev,
               pair = c("T2", "T6"),
               rank = "species",
               jitter_width_pp = 0.15,
               jitter_height_pp = 0.15)

# interactive version of the same plot
scatter_interactive(POP_prev,
               pair = c("T2", "T6"),
               rank = "species",
               jitter_width_pp = 0.15,
               jitter_height_pp = 0.15)

