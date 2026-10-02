# ================================================================
# IRN-ONLY MANTEL + PARTIAL MANTEL ANALYSIS  (simplified version)
#
# Macroinvertebrate beta diversity versus:
#   1. Environmental distance
#   2. Geographic distance
#
# Trait groups: Aerial, Purely Aquatic, Total
# Periods:      First–Sixth time
# Indices:      Jaccard (Total / JTU / JNE), Bray-Curtis (Total / BAL / GRA)
#
# Reproducibility: BASE_SEED = 20260920, N_PERM = 999
#
# Notes:
#   * IRN basin: empty samples are removed in build_comm_matrix and
#     there is NO empty-community flagging (unlike the ERN version).
#   * Ordinary and Partial Mantel share one prepared set of analysis
#     items (community, environmental, geographic, beta objects).
# ================================================================

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(openxlsx)
library(geosphere)
library(betapart)
library(ggplot2)


# ================================================================
# 1. SETTINGS
# ================================================================

file <- "FUNCAP_FWB.xlsx"

TARGET_BASIN <- "IRN"

periods <- c(
  "First time", "Second time", "Third time",
  "Fourth time", "Fifth time", "Sixth time"
)

trait_groups <- c("Aerial", "Purely Aquatic", "Total")

N_PERM <- 999
BASE_SEED <- 20260920

env_vars <- c(
  "Conductivity", "Oxygen", "pH", "Temperature",
  "TSLF", "TSLP", "TSENFlowing"
)

coord_vars <- c("lon", "lat")


# ================================================================
# 2. REPRODUCIBLE SEEDS + MANTEL WRAPPERS
# ================================================================

seed_from_labels <- function(...) {
  
  key <- paste(..., sep = "|")
  vals <- utf8ToInt(key)
  
  seed <- BASE_SEED + sum(vals * seq_along(vals))
  seed <- seed %% 2147483646
  if (seed < 1) seed <- seed + 1
  
  as.integer(seed)
}

run_mantel_reproducible <- function(
    x, y,
    method = "spearman",
    permutations = N_PERM,
    seed_key) {
  
  set.seed(seed_from_labels(seed_key))
  
  vegan::mantel(x, y, method = method, permutations = permutations)
}

run_partial_mantel_reproducible <- function(
    x, y, z,
    method = "spearman",
    permutations = N_PERM,
    seed_key) {
  
  set.seed(seed_from_labels(seed_key))
  
  vegan::mantel.partial(x, y, z, method = method,
                        permutations = permutations)
}


# ================================================================
# 3. IMPORT MACROINVERTEBRATE DATA
# ================================================================
#
# Row 1 = Class, Row 2 = Order, Row 3 = Trait, Row 4 = Family,
# Row 5 = real column names
#
# ================================================================

raw <- read_excel(
  file,
  sheet = "macroinvertebrates",
  col_names = FALSE
)

meta_cols  <- 1:6
meta_names <- c("SU", "Basin", "Flow", "Period", "Date", "Codes")

tax_class  <- unlist(raw[1, -meta_cols])
tax_order  <- unlist(raw[2, -meta_cols])
tax_traits <- unlist(raw[3, -meta_cols])
tax_family <- unlist(raw[4, -meta_cols])

dat <- raw[-c(1:4), ]

colnames(dat) <- paste0("V", seq_len(ncol(dat)))
colnames(dat)[meta_cols] <- meta_names

taxa_cols <- setdiff(colnames(dat), meta_names)

dat <- dat %>%
  mutate(across(all_of(taxa_cols),
                ~ suppressWarnings(as.numeric(.))))


# ================================================================
# 4. TAXONOMY + IRN SUBSET
# ================================================================

taxonomy_original <- data.frame(
  Taxon  = taxa_cols,
  Class  = tax_class,
  Order  = tax_order,
  Trait  = tax_traits,
  Family = tax_family,
  stringsAsFactors = FALSE
)

## standardize trait labels so "Purely_Aquatic" and "Purely Aquatic"
## are not treated as different groups
taxonomy_original <- taxonomy_original %>%
  mutate(
    Trait = trimws(as.character(Trait)),
    Trait = case_when(
      tolower(gsub("_", " ", Trait)) == "aerial" ~
        "Aerial",
      tolower(gsub("_", " ", Trait)) == "purely aquatic" ~
        "Purely Aquatic",
      TRUE ~ Trait
    )
  )

## keep only taxa occurring in the IRN
dat_irn <- dat %>%
  filter(Basin == TARGET_BASIN)

taxa_present_irn <- taxa_cols[
  colSums(dat_irn[, taxa_cols, drop = FALSE], na.rm = TRUE) > 0
]

dat_irn <- dat_irn %>%
  select(-all_of(setdiff(taxa_cols, taxa_present_irn)))

taxonomy <- taxonomy_original %>%
  filter(Taxon %in% taxa_present_irn)


# ================================================================
# 5. TRAIT GROUPS
# ================================================================

trait_taxa <- list(
  Aerial = taxonomy %>%
    filter(Trait == "Aerial") %>%
    pull(Taxon),
  
  Purely_Aquatic = taxonomy %>%
    filter(Trait == "Purely Aquatic") %>%
    pull(Taxon),
  
  Total = taxonomy %>%
    pull(Taxon)
)

for (g in trait_groups) {
  
  key <- gsub(" ", "_", g)
  
  cat(g, "taxa:", length(trait_taxa[[key]]), "\n")
  
  if (length(trait_taxa[[key]]) == 0) {
    warning(
      "No ", g, " taxa found in the IRN dataset. ",
      "Check the Trait column in the Excel file."
    )
  }
}

cat(
  "\nPurely Aquatic taxa included:\n"
)
print(trait_taxa$Purely_Aquatic)


# ================================================================
# 6. ENVIRONMENTAL DATA
# ================================================================

env_scaled <- read_excel(
  file,
  sheet = "envar"
) %>%
  filter(Basin == TARGET_BASIN)

required_env <- c("SU", "Period", env_vars, coord_vars)

missing_env <- setdiff(required_env, colnames(env_scaled))

if (length(missing_env) > 0) {
  stop(
    "Missing environmental columns: ",
    paste(missing_env, collapse = ", ")
  )
}

## scaling across the complete IRN dataset
env_scaled <- env_scaled %>%
  mutate(across(all_of(env_vars), ~ as.numeric(scale(.)))) %>%
  filter(
    !is.na(SU),
    !is.na(Period),
    complete.cases(
      across(all_of(c(env_vars, coord_vars)))
    )
  )


# ================================================================
# 7. COMMUNITY MATRIX
#
# IRN rule: completely empty samples ARE removed here; no
# empty-community flagging is needed downstream.
# ================================================================

build_comm_matrix <- function(dat_irn, taxa_names, target_period) {
  
  selected_cols <- intersect(taxa_names, colnames(dat_irn))
  
  if (length(selected_cols) == 0) return(NULL)
  
  comm <- dat_irn %>%
    filter(Period == target_period) %>%
    select(SU, Period, all_of(selected_cols)) %>%
    group_by(SU, Period) %>%
    summarise(
      across(all_of(selected_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  ## remove completely empty samples
  keep_samples <- comm %>%
    select(-SU, -Period) %>%
    rowSums(na.rm = TRUE) > 0
  
  comm <- comm[keep_samples, , drop = FALSE]
  
  if (nrow(comm) < 3) return(NULL)
  
  comm
}


# ================================================================
# 8. BETA DIVERSITY
# ================================================================

compute_beta_diversity <- function(
    mat,
    index_type = c("jaccard", "bray")) {
  
  index_type <- match.arg(index_type)
  
  if (is.null(mat)) return(NULL)
  
  abundance <- mat %>%
    select(-SU, -Period)
  
  abundance <- abundance[
    ,
    colSums(abundance) > 0,
    drop = FALSE
  ]
  
  if (ncol(abundance) == 0 || nrow(abundance) < 3) {
    return(NULL)
  }
  
  if (index_type == "jaccard") {
    
    mat_pa <- vegan::decostand(abundance, method = "pa")
    
    beta <- betapart::beta.pair(mat_pa, index.family = "jaccard")
    
    return(list(
      Jaccard_Total = beta$beta.jac,
      Jaccard_JTU   = beta$beta.jtu,
      Jaccard_JNE   = beta$beta.jne
    ))
  }
  
  beta <- betapart::beta.pair.abund(abundance, index.family = "bray")
  
  list(
    Bray_Total = beta$beta.bray,
    Bray_BAL   = beta$beta.bray.bal,
    Bray_GRA   = beta$beta.bray.gra
  )
}


# ================================================================
# 9. ENVIRONMENTAL + GEOGRAPHIC DISTANCES
# ================================================================

prepare_distance_matrices <- function(su_ids, target_period) {
  
  env_period <- env_scaled %>%
    filter(Period == target_period, SU %in% su_ids) %>%
    group_by(SU) %>%
    summarise(
      across(all_of(c(env_vars, coord_vars)), ~ first(.x)),
      .groups = "drop"
    )
  
  common_su <- intersect(su_ids, env_period$SU)
  
  if (length(common_su) < 3) return(NULL)
  
  env_period <- env_period %>%
    filter(SU %in% common_su) %>%
    arrange(match(SU, common_su))
  
  env_dist <- env_period %>%
    select(all_of(env_vars)) %>%
    as.matrix() %>%
    dist(method = "euclidean")
  
  geo_dist <- env_period %>%
    select(lon, lat) %>%
    as.matrix() %>%
    geosphere::distm(fun = geosphere::distHaversine) %>%
    as.dist()
  
  list(
    SU  = common_su,
    Env = env_dist,
    Geo = geo_dist
  )
}


# ================================================================
# 10. SHARED ANALYSIS PREPARATION
#
# Everything needed by BOTH Mantel analyses is prepared once here.
# ================================================================

prepare_analysis_items <- function() {
  
  items <- list()
  counter <- 0
  
  for (group_name in trait_groups) {
    
    taxa_names <- trait_taxa[[gsub(" ", "_", group_name)]]
    
    for (target_period in periods) {
      
      comm <- build_comm_matrix(
        dat_irn, taxa_names, target_period
      )
      
      if (is.null(comm)) next
      
      distances <- prepare_distance_matrices(
        su_ids = comm$SU,
        target_period = target_period
      )
      
      if (is.null(distances)) next
      
      common_su <- distances$SU
      
      comm <- comm %>%
        filter(SU %in% common_su) %>%
        arrange(match(SU, common_su))
      
      if (nrow(comm) < 3) next
      
      beta_list <- c(
        compute_beta_diversity(comm, "jaccard"),
        compute_beta_diversity(comm, "bray")
      )
      
      if (length(beta_list) == 0) next
      
      for (beta_index in names(beta_list)) {
        
        beta_dist <- beta_list[[beta_index]]
        
        if (is.null(beta_dist)) next
        
        if (
          attr(beta_dist, "Size") != length(common_su)
        ) {
          next
        }
        
        counter <- counter + 1
        
        items[[counter]] <- list(
          Group      = group_name,
          Period     = target_period,
          Beta_Index = beta_index,
          Beta       = beta_dist,
          Env        = distances$Env,
          Geo        = distances$Geo,
          N_Sites    = attr(beta_dist, "Size")
        )
      }
    }
  }
  
  items
}

analysis_items <- prepare_analysis_items()

cat(
  "\nPrepared", length(analysis_items),
  "group x period x index combinations\n"
)


# ================================================================
# 11. ORDINARY MANTEL
# ================================================================

run_ordinary_item <- function(item) {
  
  group_name    <- item$Group
  target_period <- item$Period
  beta_index    <- item$Beta_Index
  
  m_env <- run_mantel_reproducible(
    item$Beta,
    item$Env,
    seed_key = paste(
      "Ordinary", group_name, target_period,
      beta_index, "Env", sep = "|"
    )
  )
  
  m_geo <- run_mantel_reproducible(
    item$Beta,
    item$Geo,
    seed_key = paste(
      "Ordinary", group_name, target_period,
      beta_index, "Geo", sep = "|"
    )
  )
  
  data.frame(
    Group        = group_name,
    Period       = target_period,
    Beta_Index   = beta_index,
    N_Sites      = item$N_Sites,
    Mantel_r_Env = as.numeric(m_env$statistic),
    P_Env        = m_env$signif,
    Mantel_r_Geo = as.numeric(m_geo$statistic),
    P_Geo        = m_geo$signif,
    stringsAsFactors = FALSE
  )
}

cat("\nRunning ordinary Mantel tests...\n")

ordinary_table <- bind_rows(
  lapply(
    analysis_items,
    function(item) {
      
      message(
        "Ordinary: ", item$Group, " | ",
        item$Period, " | ", item$Beta_Index
      )
      
      run_ordinary_item(item)
    }
  )
)

print(ordinary_table)


# ================================================================
# 12. PARTIAL MANTEL
# ================================================================

run_partial_item <- function(item) {
  
  group_name    <- item$Group
  target_period <- item$Period
  beta_index    <- item$Beta_Index
  
  ## environment vs geography (one Mantel gives both r and p)
  env_geo_test <- run_mantel_reproducible(
    item$Env,
    item$Geo,
    seed_key = paste(
      "EnvGeo", group_name, target_period,
      beta_index, sep = "|"
    )
  )
  
  ## beta ~ environment | geography
  pm_env <- run_partial_mantel_reproducible(
    item$Beta,
    item$Env,
    item$Geo,
    seed_key = paste(
      "Partial", group_name, target_period,
      beta_index, "Env", sep = "|"
    )
  )
  
  ## beta ~ geography | environment
  pm_geo <- run_partial_mantel_reproducible(
    item$Beta,
    item$Geo,
    item$Env,
    seed_key = paste(
      "Partial", group_name, target_period,
      beta_index, "Geo", sep = "|"
    )
  )
  
  data.frame(
    Group         = group_name,
    Period        = target_period,
    Beta_Index    = beta_index,
    EnvGeo_r      = as.numeric(env_geo_test$statistic),
    EnvGeo_p      = env_geo_test$signif,
    Partial_r_Env = as.numeric(pm_env$statistic),
    Partial_p_Env = pm_env$signif,
    Partial_r_Geo = as.numeric(pm_geo$statistic),
    Partial_p_Geo = pm_geo$signif,
    N_Sites       = item$N_Sites,
    stringsAsFactors = FALSE
  )
}

cat("\nRunning partial Mantel tests...\n")

partial_table <- bind_rows(
  lapply(
    analysis_items,
    function(item) {
      
      message(
        "Partial: ", item$Group, " | ",
        item$Period, " | ", item$Beta_Index
      )
      
      run_partial_item(item)
    }
  )
)

print(partial_table)


# ================================================================
# 13. PROCESS CLASSIFICATION (computed ONCE, reused everywhere)
# ================================================================

partial_table <- partial_table %>%
  mutate(
    Process = case_when(
      Partial_p_Env < 0.05 & Partial_p_Geo >= 0.05 ~
        "Consistent with environmental filtering / species sorting",
      
      Partial_p_Env >= 0.05 & Partial_p_Geo < 0.05 ~
        "Consistent with dispersal limitation",
      
      Partial_p_Env < 0.05 & Partial_p_Geo < 0.05 ~
        "Joint environmental + spatial structure",
      
      TRUE ~
        "No detectable process"
    )
  )

write.xlsx(
  list(
    Ordinary_Mantel = ordinary_table,
    Partial_Mantel  = partial_table
  ),
  file = "IRN_Mantel_Results.xlsx",
  overwrite = TRUE
)


# ================================================================
# 14. PLOTTING DATA
# ================================================================

flow_status_lookup <- c(
  "First time"  = "Flowing",
  "Second time" = "Non-flowing",
  "Third time"  = "Non-flowing",
  "Fourth time" = "Flowing",
  "Fifth time"  = "Flowing",
  "Sixth time"  = "Flowing"
)

period_flow_levels <- paste0(periods, "\n(",
                             unname(flow_status_lookup[periods]), ")")

plot_data_clean <- partial_table %>%
  select(
    Group, Period, Beta_Index,
    Partial_r_Env, Partial_p_Env,
    Partial_r_Geo, Partial_p_Geo
  ) %>%
  pivot_longer(
    cols = c(
      Partial_r_Env, Partial_p_Env,
      Partial_r_Geo, Partial_p_Geo
    ),
    names_to = c(".value", "Predictor"),
    names_pattern = "Partial_(r|p)_(Env|Geo)"
  ) %>%
  mutate(
    
    Predictor = ifelse(
      Predictor == "Env",
      "Environment (geography controlled)",
      "Geography (environment controlled)"
    ),
    
    Signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    
    Group = factor(
      Group,
      levels = c("Aerial", "Purely Aquatic", "Total")
    ),
    
    Period = factor(Period, levels = periods),
    
    Period_Flow = factor(
      paste0(
        Period,
        "\n(",
        unname(flow_status_lookup[as.character(Period)]),
        ")"
      ),
      levels = period_flow_levels
    ),
    
    Beta_Label = factor(
      dplyr::recode(
        Beta_Index,
        "Bray_Total"    = "B-Tot",
        "Bray_BAL"      = "BAL",
        "Bray_GRA"      = "GRA",
        "Jaccard_Total" = "J-Tot",
        "Jaccard_JTU"   = "JTU",
        "Jaccard_JNE"   = "JNE"
      ),
      levels = c("B-Tot", "BAL", "GRA", "J-Tot", "JTU", "JNE")
    ),
    
    Predictor = factor(
      Predictor,
      levels = c(
        "Environment (geography controlled)",
        "Geography (environment controlled)"
      )
    ),
    
    ## significance-label positions
    y_label = ifelse(r >= 0, r + 0.04, r - 0.04),
    
    vjust_lbl = ifelse(r >= 0, 0, 1)
  ) %>%
  filter(!is.na(Beta_Label), !is.na(r))


# ================================================================
# 15. PARTIAL MANTEL FIGURE
# ================================================================

p4_base1 <- ggplot(
  plot_data_clean,
  aes(x = Beta_Label, y = r, fill = Predictor)
) +
  
  geom_col(
    position = position_dodge(width = 0.8, preserve = "single"),
    width = 0.75,
    color = "grey20",
    linewidth = 0.3
  ) +
  
  geom_text(
    aes(label = Signif, y = y_label, vjust = vjust_lbl),
    position = position_dodge(width = 0.8, preserve = "single"),
    size = 2.6,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey20") +
  
  facet_grid(
    Group ~ Period_Flow,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Environment (geography controlled)" = "#4D4D4D",
      "Geography (environment controlled)" = "#C8C8C8"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  
  labs(
    x = NULL,
    y = "Partial Mantel correlation (r)",
    fill = NULL
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(
      face = "bold", size = 9, color = "white"
    ),
    strip.background = element_rect(
      fill = "grey25", color = NA
    ),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.3,
      linetype = "dashed"
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.5)
  )

print(p4_base1)

ggsave(
  filename = "IRN_Partial_Mantel.tiff",
  plot = p4_base1,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ================================================================
# 16. SUMMARY + SIGNIFICANT-ONLY EXPORTS
# ================================================================

summary_table <- partial_table %>%
  select(
    Group,
    Period,
    Beta_Index,
    N_Sites,
    Partial_r_Env,
    Partial_p_Env,
    Partial_r_Geo,
    Partial_p_Geo,
    Process
  )

print(summary_table)

write.xlsx(
  summary_table,
  file = "IRN_Partial_Mantel_Summary.xlsx",
  overwrite = TRUE
)

significant_partial <- partial_table %>%
  filter(
    Partial_p_Env < 0.05 |
      Partial_p_Geo < 0.05
  ) %>%
  arrange(Group, Period, Beta_Index)

print(significant_partial)

write.xlsx(
  significant_partial,
  file = "IRN_Partial_Mantel_Significant.xlsx",
  overwrite = TRUE
)


# ================================================================
# 17. FINAL CHECKS
# ================================================================

cat(
  "\n================================================\n",
  "IRN Mantel analysis completed.\n",
  "Prepared analysis items:", length(analysis_items), "\n",
  "Ordinary Mantel results:", nrow(ordinary_table), "rows\n",
  "Partial Mantel results:", nrow(partial_table), "rows\n",
  "Significant partial Mantel results:",
  nrow(significant_partial), "rows\n",
  "================================================\n"
)

# ================================================================
# END
# ================================================================
