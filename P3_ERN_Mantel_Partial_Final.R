# ================================================================
##  P2: Testing prediction P3 from FWB manuscript (Menezes et al.)
## Figure: Figure 6 
## ERN-ONLY MANTEL + PARTIAL MANTEL ANALYSIS  (simplified version)
#
# Macroinvertebrate beta diversity versus:
#   1. Environmental distance
#   2. Geographic distance
#
# Trait groups: Aerial, Purely Aquatic, Total Macroinvertebrates
# Periods:      First, Fifth, Sixth time
# Indices:      Jaccard (Total / JTU / JNE), Bray-Curtis (Total / BAL / GRA)
#
# Analyses:
#   Ordinary Mantel:  beta ~ environment ;  beta ~ geography
#   Partial Mantel:   beta ~ environment | geography
#                     beta ~ geography    | environment
#
# COMMUNITY-MATRIX RULE:
#   Total Jaccard / Bray-Curtis: empty vs non-empty = 1, empty vs empty = 0
#   Components (JTU/JNE, BAL/GRA): NA for comparisons involving empty
#   communities -> tests flagged "Not estimable: empty community" (shown
#   as † in the figure).
#
# Reproducibility: BASE_SEED = 20260920, N_PERM = 999
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

TARGET_BASIN <- "ERN"

periods <- c("First time", "Fifth time", "Sixth time")

trait_groups <- c("Aerial", "Purely_Aquatic", "Total")

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
# 4. TAXONOMY TABLE + ERN SUBSET
# ================================================================

taxonomy_original <- data.frame(
  Taxon = taxa_cols,
  Class = tax_class,
  Order = tax_order,
  Trait = tax_traits,
  Family = tax_family,
  stringsAsFactors = FALSE
)

dat_ern <- dat %>%
  filter(Basin == TARGET_BASIN)

## drop taxa absent from the ERN
taxa_present_ern <- taxa_cols[
  colSums(dat_ern[, taxa_cols, drop = FALSE], na.rm = TRUE) > 0
]

dat_ern <- dat_ern %>%
  select(-all_of(setdiff(taxa_cols, taxa_present_ern)))

taxonomy <- taxonomy_original %>%
  filter(Taxon %in% taxa_present_ern)


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
  
  cat(
    g, "taxa:",
    length(trait_taxa[[g]]),
    "\n"
  )
  
  if (length(trait_taxa[[g]]) == 0) {
    stop("No taxa found for group ", g,
         " - check the Trait column.")
  }
}


# ================================================================
# 6. ENVIRONMENTAL DATA
# ================================================================

env_scaled <- read_excel(
  file,
  sheet = "envar"
) %>%
  filter(Basin == TARGET_BASIN) %>%
  mutate(across(all_of(env_vars), ~ as.numeric(scale(.)))) %>%
  filter(
    !is.na(SU),
    !is.na(Period)
  ) %>%
  filter(
    complete.cases(
      across(all_of(c(env_vars, coord_vars)))
    )
  )

required_env <- c("SU", "Period", env_vars, coord_vars)

missing_env <- setdiff(required_env, colnames(env_scaled))

if (length(missing_env) > 0) {
  stop(
    "Missing environmental columns: ",
    paste(missing_env, collapse = ", ")
  )
}


# ================================================================
# 7. COMMUNITY MATRIX FUNCTION
#
# Completely empty samples are NOT removed: a sampled SU with zero
# abundance for the focal trait group is a valid sampled site and must
# remain represented in the environmental and geographic matrices.
# ================================================================

build_comm_matrix <- function(dat_ern, taxa_names, target_period) {
  
  selected_cols <- intersect(taxa_names, colnames(dat_ern))
  
  if (length(selected_cols) == 0) {
    return(NULL)
  }
  
  comm <- dat_ern %>%
    filter(Period == target_period) %>%
    select(SU, Period, all_of(selected_cols)) %>%
    group_by(SU, Period) %>%
    summarise(
      across(all_of(selected_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  if (nrow(comm) < 2) {
    return(NULL)
  }
  
  comm %>%
    mutate(
      Empty_Sample =
        rowSums(
          select(., -SU, -Period),
          na.rm = TRUE
        ) == 0
    )
}


# ================================================================
# 8. BETA-DIVERSITY FUNCTION
#
# Total indices: computed for ALL sampled SUs, including empty samples
#   (empty vs non-empty = 1, empty vs empty = 0)
# Component indices (JTU/JNE, BAL/GRA): only among non-empty samples;
#   comparisons involving empty samples are NA.
# ================================================================

compute_beta_diversity <- function(mat, index_type = "jaccard") {
  
  if (is.null(mat)) {
    return(NULL)
  }
  
  mat_abund <- mat %>%
    select(-SU, -Period, -Empty_Sample)
  
  ## remove taxa absent from all retained samples
  mat_abund <- mat_abund[
    ,
    colSums(mat_abund, na.rm = TRUE) > 0,
    drop = FALSE
  ]
  
  if (nrow(mat_abund) < 2) {
    return(NULL)
  }
  
  empty_samples <- rowSums(mat_abund, na.rm = TRUE) == 0
  n <- nrow(mat_abund)
  
  
  ## ----------------------------------------------------------
  ## shared tail: component matrices on non-empty samples only
  ## ----------------------------------------------------------
  
  compute_components <- function(beta_nonempty) {
    
    comp1 <- matrix(NA_real_, nrow = n, ncol = n)
    comp2 <- matrix(NA_real_, nrow = n, ncol = n)
    
    diag(comp1) <- 0
    diag(comp2) <- 0
    
    nonempty <- which(!empty_samples)
    
    if (length(nonempty) >= 2) {
      
      m1 <- as.matrix(beta_nonempty[[1]])
      m2 <- as.matrix(beta_nonempty[[2]])
      
      comp1[nonempty, nonempty] <- m1
      comp2[nonempty, nonempty] <- m2
    }
    
    list(as.dist(comp1), as.dist(comp2))
  }
  
  
  ## ----------------------------------------------------------
  ## JACCARD
  ## ----------------------------------------------------------
  
  if (index_type == "jaccard") {
    
    mat_pa <- vegan::decostand(mat_abund, method = "pa")
    
    ## total Jaccard (manual loop: handles empty samples)
    jac_total <- matrix(0, nrow = n, ncol = n)
    
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        
        if (empty_samples[i] && empty_samples[j]) {
          
          d <- 0
          
        } else if (empty_samples[i] || empty_samples[j]) {
          
          d <- 1
          
        } else {
          
          a <- sum(mat_pa[i, ] == 1 & mat_pa[j, ] == 1)
          b <- sum(mat_pa[i, ] == 1 & mat_pa[j, ] == 0)
          c <- sum(mat_pa[i, ] == 0 & mat_pa[j, ] == 1)
          
          denominator <- a + b + c
          
          d <- ifelse(
            denominator == 0,
            0,
            (b + c) / denominator
          )
        }
        
        jac_total[i, j] <- d
        jac_total[j, i] <- d
      }
    }
    
    nonempty <- which(!empty_samples)
    
    beta_nonempty <- NULL
    
    if (length(nonempty) >= 2) {
      
      beta_nonempty <- betapart::beta.pair(
        mat_pa[nonempty, , drop = FALSE],
        index.family = "jaccard"
      )
    }
    
    comps <- list(NULL, NULL)
    
    if (!is.null(beta_nonempty)) {
      
      comps <- compute_components(
        list(
          beta_nonempty$beta.jtu,
          beta_nonempty$beta.jne
        )
      )
    }
    
    return(
      list(
        Total        = as.dist(jac_total),
        JTU          = comps[[1]],
        JNE          = comps[[2]],
        Empty_Samples = empty_samples
      )
    )
  }
  
  
  ## ----------------------------------------------------------
  ## BRAY-CURTIS
  ## ----------------------------------------------------------
  
  if (index_type == "bray") {
    
    bray_total <- matrix(0, nrow = n, ncol = n)
    
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        
        xi <- mat_abund[i, ]
        xj <- mat_abund[j, ]
        
        si <- sum(xi)
        sj <- sum(xj)
        
        if (si == 0 && sj == 0) {
          
          d <- 0
          
        } else if (si == 0 || sj == 0) {
          
          d <- 1
          
        } else {
          
          d <- sum(abs(xi - xj)) / sum(xi + xj)
        }
        
        bray_total[i, j] <- d
        bray_total[j, i] <- d
      }
    }
    
    nonempty <- which(!empty_samples)
    
    beta_nonempty <- NULL
    
    if (length(nonempty) >= 2) {
      
      beta_nonempty <- betapart::beta.pair.abund(
        mat_abund[nonempty, , drop = FALSE],
        index.family = "bray"
      )
    }
    
    comps <- list(NULL, NULL)
    
    if (!is.null(beta_nonempty)) {
      
      comps <- compute_components(
        list(
          beta_nonempty$beta.bray.bal,
          beta_nonempty$beta.bray.gra
        )
      )
    }
    
    return(
      list(
        Total         = as.dist(bray_total),
        BAL           = comps[[1]],
        GRA           = comps[[2]],
        Empty_Samples = empty_samples
      )
    )
  }
  
  NULL
}


# ================================================================
# 9. ENVIRONMENTAL + GEOGRAPHIC DISTANCE MATRICES
# ================================================================

prepare_distance_matrices <- function(su_ids, target_period) {
  
  env_period <- env_scaled %>%
    filter(Period == target_period, SU %in% su_ids) %>%
    group_by(SU) %>%
    summarise(
      across(
        all_of(c(env_vars, coord_vars)),
        ~ first(.x)
      ),
      .groups = "drop"
    )
  
  ## match SUs between community and environmental data
  common_su <- intersect(su_ids, env_period$SU)
  
  if (length(common_su) < 3) {
    return(NULL)
  }
  
  env_period <- env_period %>%
    filter(SU %in% common_su) %>%
    arrange(match(SU, common_su))
  
  ## environmental distance
  env_dist <- env_period %>%
    select(all_of(env_vars)) %>%
    as.matrix() %>%
    dist(method = "euclidean")
  
  ## geographic distance (haversine, metres)
  geo_dist_m <- env_period %>%
    select(lon, lat) %>%
    as.matrix() %>%
    geosphere::distm(fun = geosphere::distHaversine)
  
  list(
    SU  = common_su,
    Env = env_dist,
    Geo = as.dist(geo_dist_m)
  )
}


# ================================================================
# 10. SHARED PREPARATION FOR BOTH ANALYSES
#
# Returns a flat list of per-(group, period, index) items holding the
# beta distance matrix, env/geo distances, status and N_empty — or NULL
# when the combination cannot be analysed.
# ================================================================

prepare_analysis_items <- function() {
  
  items <- list()
  
  for (group_name in trait_groups) {
    
    taxa_names <- trait_taxa[[group_name]]
    
    for (target_period in periods) {
      
      comm <- build_comm_matrix(
        dat_ern, taxa_names, target_period
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
      
      n_empty <- sum(comm$Empty_Sample)
      
      beta_jaccard <- compute_beta_diversity(comm, "jaccard")
      beta_bray    <- compute_beta_diversity(comm, "bray")
      
      
      for (beta_list in list(beta_jaccard, beta_bray)) {
        
        if (is.null(beta_list)) next
        
        index_prefix <- ifelse(
          identical(beta_list, beta_jaccard),
          "Jaccard_",
          "Bray_"
        )
        
        for (idx_name in names(beta_list)) {
          
          if (idx_name == "Empty_Samples") next
          
          beta_dist <- beta_list[[idx_name]]
          
          if (
            attr(beta_dist, "Size") != length(common_su)
          ) {
            next
          }
          
          beta_values <- as.vector(beta_dist)
          
          status <- ifelse(
            any(!is.finite(beta_values)),
            "Not estimable: empty community",
            "Estimated"
          )
          
          items[[length(items) + 1]] <- list(
            group      = group_name,
            period     = target_period,
            index_name = paste0(index_prefix, idx_name),
            beta_dist  = beta_dist,
            env_dist   = distances$Env,
            geo_dist   = distances$Geo,
            status     = status,
            n_empty    = n_empty
          )
        }
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

run_ordinary_result <- function(item) {
  
  if (item$status != "Estimated") {
    
    return(
      data.frame(
        Group      = item$group,
        Period     = item$period,
        Beta_Index = item$index_name,
        N_Sites    = attr(item$beta_dist, "Size"),
        N_Empty    = item$n_empty,
        Mantel_r_Env = NA_real_,
        P_Env        = NA_real_,
        Mantel_r_Geo = NA_real_,
        P_Geo        = NA_real_,
        Status       = item$status,
        stringsAsFactors = FALSE
      )
    )
  }
  
  m_env <- run_mantel_reproducible(
    item$beta_dist,
    item$env_dist,
    seed_key = paste(
      "Ordinary", item$group, item$period,
      item$index_name, "Env", sep = "|"
    )
  )
  
  m_geo <- run_mantel_reproducible(
    item$beta_dist,
    item$geo_dist,
    seed_key = paste(
      "Ordinary", item$group, item$period,
      item$index_name, "Geo", sep = "|"
    )
  )
  
  data.frame(
    Group      = item$group,
    Period     = item$period,
    Beta_Index = item$index_name,
    N_Sites    = attr(item$beta_dist, "Size"),
    N_Empty    = item$n_empty,
    Mantel_r_Env = as.numeric(m_env$statistic),
    P_Env        = m_env$signif,
    Mantel_r_Geo = as.numeric(m_geo$statistic),
    P_Geo        = m_geo$signif,
    Status       = "Estimated",
    stringsAsFactors = FALSE
  )
}

cat("\nRunning ordinary Mantel tests...\n")

ordinary_table <- bind_rows(
  lapply(
    analysis_items,
    function(item) {
      
      message(
        "Ordinary: ", item$group, " | ",
        item$period, " | ", item$index_name
      )
      
      run_ordinary_result(item)
    }
  )
)


# ================================================================
# 12. PARTIAL MANTEL
# ================================================================

run_partial_result <- function(item) {
  
  if (item$status != "Estimated") {
    
    return(
      data.frame(
        Group      = item$group,
        Period     = item$period,
        Beta_Index = item$index_name,
        EnvGeo_r   = NA_real_,
        EnvGeo_p   = NA_real_,
        Partial_r_Env = NA_real_,
        Partial_p_Env = NA_real_,
        Partial_r_Geo = NA_real_,
        Partial_p_Geo = NA_real_,
        N_Sites    = attr(item$beta_dist, "Size"),
        N_Empty    = item$n_empty,
        Status     = item$status,
        stringsAsFactors = FALSE
      )
    )
  }
  
  ## environment vs geography
  pm_env <- run_partial_mantel_reproducible(
    item$beta_dist,
    item$env_dist,
    item$geo_dist,
    seed_key = paste(
      "Partial", item$group, item$period,
      item$index_name, "Env", sep = "|"
    )
  )
  
  ## geography vs environment
  pm_geo <- run_partial_mantel_reproducible(
    item$beta_dist,
    item$geo_dist,
    item$env_dist,
    seed_key = paste(
      "Partial", item$group, item$period,
      item$index_name, "Geo", sep = "|"
    )
  )
  
  ## env vs geo correlation (one Mantel test gives both r and p)
  env_geo <- run_mantel_reproducible(
    item$env_dist,
    item$geo_dist,
    seed_key = paste(
      "EnvGeo", item$group, item$period,
      item$index_name, sep = "|"
    )
  )
  
  data.frame(
    Group      = item$group,
    Period     = item$period,
    Beta_Index = item$index_name,
    EnvGeo_r   = as.numeric(env_geo$statistic),
    EnvGeo_p   = env_geo$signif,
    Partial_r_Env = as.numeric(pm_env$statistic),
    Partial_p_Env = pm_env$signif,
    Partial_r_Geo = as.numeric(pm_geo$statistic),
    Partial_p_Geo = pm_geo$signif,
    N_Sites    = attr(item$beta_dist, "Size"),
    N_Empty    = item$n_empty,
    Status     = "Estimated",
    stringsAsFactors = FALSE
  )
}

cat("\nRunning partial Mantel tests...\n")

partial_table <- bind_rows(
  lapply(
    analysis_items,
    function(item) {
      
      message(
        "Partial: ", item$group, " | ",
        item$period, " | ", item$index_name
      )
      
      run_partial_result(item)
    }
  )
)

print(partial_table)


# ================================================================
# 13. PROCESS CLASSIFICATION + RAW EXPORTS
# ================================================================

partial_table <- partial_table %>%
  mutate(
    Env_Significant = !is.na(Partial_p_Env) &
      Partial_p_Env < 0.05,
    Geo_Significant = !is.na(Partial_p_Geo) &
      Partial_p_Geo < 0.05,
    
    Process = case_when(
      Status != "Estimated" ~
        "Not estimable",
      Env_Significant & !Geo_Significant ~
        "Consistent with environmental filtering / species sorting",
      !Env_Significant & Geo_Significant ~
        "Consistent with dispersal limitation",
      Env_Significant & Geo_Significant ~
        "Joint environmental + spatial structure",
      TRUE ~
        "No detectable environmental or spatial structure"
    )
  )

write.xlsx(
  list(
    Ordinary_Mantel = ordinary_table,
    Partial_Mantel  = partial_table
  ),
  file = "ERN_Mantel_Results.xlsx",
  overwrite = TRUE
)


# ================================================================
# 14. PLOTTING DATA
# ================================================================

plot_data_clean <- partial_table %>%
  select(
    Group, Period, Beta_Index,
    Partial_r_Env, Partial_p_Env,
    Partial_r_Geo, Partial_p_Geo,
    Status
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
    
    ## † = not estimable
    Signif = case_when(
      Status != "Estimated" ~ "†",
      is.na(p)              ~ "",
      p < 0.001             ~ "***",
      p < 0.01              ~ "**",
      p < 0.05              ~ "*",
      TRUE                  ~ ""
    ),
    
    Group = ifelse(
      Group == "Purely_Aquatic",
      "Purely Aquatic",
      Group
    ),
    
    Group = factor(
      Group,
      levels = c("Aerial", "Purely Aquatic", "Total")
    ),
    
    Period = factor(
      Period,
      levels = periods
    ),
    
    Period_Flow = factor(
      paste0(Period, "\n(Non-flowing)"),
      levels = paste0(periods, "\n(Non-flowing)")
    ),
    
    Beta_Label = dplyr::recode(
      Beta_Index,
      "Bray_Total"    = "B-Tot",
      "Bray_BAL"      = "BAL",
      "Bray_GRA"      = "GRA",
      "Jaccard_Total" = "J-Tot",
      "Jaccard_JTU"   = "JTU",
      "Jaccard_JNE"   = "JNE"
    ),
    
    Beta_Label = factor(
      Beta_Label,
      levels = c(
        "B-Tot", "BAL", "GRA",
        "J-Tot", "JTU", "JNE"
      )
    ),
    
    Predictor = factor(
      Predictor,
      levels = c(
        "Environment (geography controlled)",
        "Geography (environment controlled)"
      )
    ),
    
    ## significance-label positions
    y_label = case_when(
      Status != "Estimated" ~ 0.08,
      r >= 0                ~ r + 0.04,
      TRUE                  ~ r - 0.04
    ),
    
    vjust_lbl = case_when(
      Status != "Estimated" ~ 0,
      r >= 0                ~ 0,
      TRUE                  ~ 1
    )
  ) %>%
  filter(!is.na(Beta_Label))

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
    linewidth = 0.3,
    na.rm = TRUE
  ) +
  
  geom_text(
    aes(label = Signif, y = y_label, vjust = vjust_lbl),
    position = position_dodge(width = 0.8, preserve = "single"),
    size = 2.8,
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
  filename = "ERN_Partial_Mantel.tiff",
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
  mutate(
    Group = ifelse(
      Group == "Purely_Aquatic",
      "Purely Aquatic",
      Group
    )
  ) %>%
  select(
    Group, Period, Beta_Index, N_Sites, N_Empty,
    EnvGeo_r, EnvGeo_p,
    Partial_r_Env, Partial_p_Env,
    Partial_r_Geo, Partial_p_Geo,
    Status, Process
  )

print(summary_table)

write.xlsx(
  summary_table,
  file = "ERN_Partial_Mantel_Summary.xlsx",
  overwrite = TRUE
)

significant_partial <- partial_table %>%
  mutate(
    Group = ifelse(
      Group == "Purely_Aquatic",
      "Purely Aquatic",
      Group
    )
  ) %>%
  filter(
    Status == "Estimated",
    Partial_p_Env < 0.05 |
      Partial_p_Geo < 0.05
  ) %>%
  arrange(Group, Period, Beta_Index)

print(significant_partial)

write.xlsx(
  significant_partial,
  file = "ERN_Partial_Mantel_Significant.xlsx",
  overwrite = TRUE
)


# ================================================================
# 17. EMPTY-SAMPLE DIAGNOSTIC
# ================================================================

empty_sample_table <- bind_rows(
  lapply(
    trait_groups,
    function(group_name) {
      
      bind_rows(
        lapply(
          periods,
          function(target_period) {
            
            comm_check <- build_comm_matrix(
              dat_ern,
              trait_taxa[[group_name]],
              target_period
            )
            
            if (is.null(comm_check)) {
              return(NULL)
            }
            
            empty_sus <- comm_check %>%
              filter(Empty_Sample) %>%
              pull(SU)
            
            data.frame(
              Group = group_name,
              Period = target_period,
              N_Sites = nrow(comm_check),
              N_Empty = length(empty_sus),
              Empty_SUs = ifelse(
                length(empty_sus) == 0,
                "",
                paste(empty_sus, collapse = ", ")
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

cat("\n================================================\n")
cat("Empty-community diagnostic:\n")
print(empty_sample_table)

write.xlsx(
  empty_sample_table,
  file = "ERN_Mantel_Empty_Sample_Diagnostic.xlsx",
  overwrite = TRUE
)


# ================================================================
# 18. FINAL CHECKS
# ================================================================

cat("\n================================================\n")
cat("ERN Mantel analysis completed.\n")
cat("Ordinary Mantel results:",  nrow(ordinary_table),    "rows\n")
cat("Partial Mantel results:",   nrow(partial_table),     "rows\n")
cat("Significant partial results:", nrow(significant_partial), "rows\n")
cat(
  "Not-estimable partial results:",
  sum(partial_table$Status != "Estimated"),
  "rows\n"
)
cat("================================================\n")

# ================================================================
# END
# ================================================================
