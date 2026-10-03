## MENEZES et al. (Freshwater Biology manuscript)

################################################################################
##
##  DESCRIPTION:
##  Codes for testing P2 and generating FIGURE 4, that spatial β-diversity is 
##  higher in ERN than in IRN whereas temporal β-diversity is lower in ERN than
##  in IRN
##  
################################################################################
## 
##  SPATIAL BETA DIVERSITY: ERN vs IRN
##  Site-level permutation test 
##
##  Permutation scheme:
##    1. Randomizes Basin labels among sites
##    2. Keeps exactly 9 sites in pseudo-ERN, 18 in pseudo-IRN
##    3. Keeps all sampling occasions of a site together
##    4. Recalculates spatial beta within each occasion, averages
##    5. Builds a null distribution of ERN - IRN differences
##
##  NOTE: no multiple-test correction is applied; reported p-values
##  are raw two-sided permutation p-values.
##
################################################################################

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(betapart) # provides functions for partitioning beta diversity into its main components, particularly species turnover and nestedness-resultant dissimilarity, based on presence–absence or abundance data. It is used here for the analyses of differences in community composition among sites.


################################################################################
## 1. READ DATA
################################################################################

file <- "FUNCAP_FWB.xlsx"

df_raw <- read_excel(
  file,
  sheet = "macroinvertebrates",
  col_names = FALSE
)

## Header structure:
##   Row 1 = Class, Row 2 = Order,
##   Row 3 = Functional group, Row 4 = column names

colnames(df_raw) <- as.character(df_raw[4, ])
df <- df_raw[-c(1:4), ]

meta_cols <- c("SU", "Basin", "Flow", "Period", "Date", "Codes")
taxa_cols <- colnames(df)[!colnames(df) %in% meta_cols]

fg_df <- data.frame(
  taxa  = taxa_cols,
  group = as.character(df_raw[3, match(taxa_cols, colnames(df_raw))]),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(group))

aerial_taxa  <- fg_df$taxa[fg_df$group == "Aerial"]
aquatic_taxa <- fg_df$taxa[fg_df$group == "Purely Aquatic"]

df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)


################################################################################
## 2. COMMUNITY DATASETS
################################################################################

## Aquatic excludes two zero-abundance SUs; existing sampling
## structure is retained (no requirement of complete site x period grids)

macro_total <- df[, c(meta_cols, taxa_cols)] %>% drop_na()
macro_aerial <- df[, c(meta_cols, aerial_taxa)] %>% drop_na()
macro_aquatic <- df[, c(meta_cols, aquatic_taxa)] %>%
  filter(!SU %in% c("BN07-C5", "BN10-C5")) %>%
  drop_na()

comm_total   <- macro_total[, taxa_cols, drop = FALSE]
comm_aerial  <- macro_aerial[, aerial_taxa, drop = FALSE]
comm_aquatic <- macro_aquatic[, aquatic_taxa, drop = FALSE]

add_site_id <- function(m) {
  m %>%
    mutate(Site_ID = sub("-C[0-9]+$", "", as.character(SU)))
}

metadata_total   <- macro_total[, meta_cols, drop = FALSE]   %>% add_site_id()
metadata_aerial  <- macro_aerial[, meta_cols, drop = FALSE]  %>% add_site_id()
metadata_aquatic <- macro_aquatic[, meta_cols, drop = FALSE] %>% add_site_id()

pa_total   <- decostand(comm_total,   method = "pa")
pa_aerial  <- decostand(comm_aerial,  method = "pa")
pa_aquatic <- decostand(comm_aquatic, method = "pa")

cat("\nSites per basin:\n")
print(metadata_total %>% distinct(Basin, Site_ID) %>% count(Basin))


################################################################################
## 3. CLEAN COMMUNITY
################################################################################

clean_community <- function(comm_matrix, metadata) {
  
  keep_rows <- rowSums(comm_matrix, na.rm = TRUE) > 0
  comm_clean <- comm_matrix[keep_rows, , drop = FALSE]
  meta_clean <- metadata[keep_rows, , drop = FALSE]
  
  keep_cols <- colSums(comm_clean, na.rm = TRUE) > 0
  comm_clean <- comm_clean[, keep_cols, drop = FALSE]
  
  if (any(is.na(comm_clean))) {
    warning("NA values found - replaced with 0")
    comm_clean[is.na(comm_clean)] <- 0
  }
  
  if (nrow(comm_clean) < 2) stop("Fewer than 2 sampling units remain!")
  if (ncol(comm_clean) < 2) stop("Fewer than 2 taxa remain!")
  
  list(comm = comm_clean, meta = meta_clean)
}


################################################################################
## 4. SPATIAL BETA WITHIN ONE SAMPLING OCCASION
##
## Beta is calculated only among DIFFERENT sites sampled during the
## SAME occasion (within-site temporal pairs are excluded).
################################################################################

calculate_spatial_beta <- function(comm_matrix, metadata,
                                   basin_assignment,
                                   method = "jaccard",
                                   component = "Total") {
  
  metadata_work <- metadata %>%
    mutate(Perm_Basin = basin_assignment[Site_ID])
  
  periods <- unique(metadata_work$Period)
  periods <- periods[!is.na(periods)]
  
  period_results <- list()
  
  for (p in periods) {
    
    idx_period <- which(metadata_work$Period == p)
    comm_period <- comm_matrix[idx_period, , drop = FALSE]
    meta_period <- metadata_work[idx_period, , drop = FALSE]
    
    keep_rows <- rowSums(comm_period, na.rm = TRUE) > 0
    comm_period <- comm_period[keep_rows, , drop = FALSE]
    meta_period <- meta_period[keep_rows, , drop = FALSE]
    
    for (basin in c("ERN", "IRN")) {
      
      idx_basin <- which(meta_period$Perm_Basin == basin)
      comm_basin <- comm_period[idx_basin, , drop = FALSE]
      
      if (nrow(comm_basin) < 2) {
        period_results[[length(period_results) + 1]] <-
          data.frame(Period = p, Basin = basin,
                     Beta = NA_real_, SD = NA_real_,
                     N_Pairs = NA_integer_, SE = NA_real_)
        next
      }
      
      keep_species <- colSums(comm_basin, na.rm = TRUE) > 0
      comm_basin <- comm_basin[, keep_species, drop = FALSE]
      
      if (ncol(comm_basin) < 2) {
        period_results[[length(period_results) + 1]] <-
          data.frame(Period = p, Basin = basin,
                     Beta = NA_real_, SD = NA_real_,
                     N_Pairs = NA_integer_, SE = NA_real_)
        next
      }
      
      if (method == "jaccard") {
        
        comm_basin <- decostand(comm_basin, method = "pa")
        
        beta <- beta.pair(comm_basin, index.family = "jaccard")
        
        beta_values <- switch(
          component,
          Total      = beta$beta.jac,
          Turnover   = beta$beta.jtu,
          Nestedness = beta$beta.jne
        )
      }
      
      if (method == "bray") {
        
        beta <- beta.pair.abund(comm_basin, index.family = "bray")
        
        beta_values <- switch(
          component,
          Total      = beta$beta.bray,
          Turnover   = beta$beta.bray.bal,
          Nestedness = beta$beta.bray.gra
        )
      }
      
      beta_numeric <- as.numeric(beta_values)
      
      mean_beta <- mean(beta_numeric, na.rm = TRUE)
      sd_beta   <- sd(beta_numeric, na.rm = TRUE)
      n_pairs   <- sum(!is.na(beta_numeric))
      
      period_results[[length(period_results) + 1]] <-
        data.frame(
          Period  = p,
          Basin   = basin,
          Beta    = mean_beta,
          SD      = sd_beta,
          N_Pairs = n_pairs,
          SE      = sd_beta / sqrt(n_pairs)
        )
    }
  }
  
  period_results <- bind_rows(period_results)
  
  overall_results <- period_results %>%
    group_by(Basin) %>%
    summarise(
      Mean_Beta = mean(Beta, na.rm = TRUE),
      N_Periods = sum(!is.na(Beta)),
      .groups = "drop"
    )
  
  list(period = period_results, overall = overall_results)
}


################################################################################
## 5. SITE-LEVEL PERMUTATION TEST
################################################################################

permutation_spatial_beta <- function(comm_matrix, metadata,
                                     method = "jaccard",
                                     component = "Total",
                                     nperm = 4999,
                                     seed = 123) {
  
  set.seed(seed)
  
  if (!"Site_ID" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Site_ID = sub("-C[0-9]+$", "", as.character(SU)))
  }
  
  cleaned <- clean_community(comm_matrix, metadata)
  comm_clean <- cleaned$comm
  meta_clean <- cleaned$meta
  
  site_info <- meta_clean %>%
    distinct(Site_ID, Basin)
  
  if (!all(c("ERN", "IRN") %in% unique(site_info$Basin))) {
    stop("Both ERN and IRN must be present.")
  }
  
  n_ERN <- sum(site_info$Basin == "ERN")
  n_IRN <- sum(site_info$Basin == "IRN")
  
  observed_assignment <- site_info$Basin
  names(observed_assignment) <- site_info$Site_ID
  
  observed <- calculate_spatial_beta(
    comm_matrix = comm_clean,
    metadata = meta_clean,
    basin_assignment = observed_assignment,
    method = method,
    component = component
  )
  
  observed_values <- observed$overall
  
  observed_ERN <- observed_values$Mean_Beta[
    observed_values$Basin == "ERN"]
  observed_IRN <- observed_values$Mean_Beta[
    observed_values$Basin == "IRN"]
  
  observed_difference <- observed_ERN - observed_IRN
  
  null_difference <- numeric(nperm)
  site_ids <- site_info$Site_ID
  
  cat("\nRunning", nperm, "permutations:",
      method, "-", component, "\n")
  
  for (i in seq_len(nperm)) {
    
    ## pseudo-ERN = random subset of exactly n_ERN sites
    pseudo_ERN <- sample(site_ids, size = n_ERN, replace = FALSE)
    
    perm_assignment <- rep("IRN", length(site_ids))
    names(perm_assignment) <- site_ids
    perm_assignment[pseudo_ERN] <- "ERN"
    
    perm_result <- calculate_spatial_beta(
      comm_matrix = comm_clean,
      metadata = meta_clean,
      basin_assignment = perm_assignment,
      method = method,
      component = component
    )
    
    perm_values <- perm_result$overall
    
    perm_ERN <- perm_values$Mean_Beta[perm_values$Basin == "ERN"]
    perm_IRN <- perm_values$Mean_Beta[perm_values$Basin == "IRN"]
    
    null_difference[i] <- perm_ERN - perm_IRN
  }
  
  ## two-sided permutation p-value (with +1 correction)
  p_value <- (
    sum(abs(null_difference) >= abs(observed_difference), na.rm = TRUE) + 1
  ) / (nperm + 1)
  
  ## one-sided p-value (ERN > IRN)
  p_value_ERN_greater <- (
    sum(null_difference >= observed_difference, na.rm = TRUE) + 1
  ) / (nperm + 1)
  
  list(
    method              = method,
    component           = component,
    nperm               = nperm,
    n_ERN               = n_ERN,
    n_IRN               = n_IRN,
    observed_ERN        = observed_ERN,
    observed_IRN        = observed_IRN,
    observed_difference = observed_difference,
    p_value_two_sided   = p_value,
    p_value_ERN_greater = p_value_ERN_greater,
    null_difference     = null_difference,
    observed_period     = observed$period
  )
}


################################################################################
## 6. RUN ALL 18 TESTS
##    3 communities x 2 indices x 3 components
################################################################################

run_all_spatial_permutations <- function(nperm = 4999, seed = 123) {
  
  community_data <- list(
    Total = list(
      comm_jaccard = pa_total,
      comm_bray    = comm_total,
      metadata     = metadata_total
    ),
    Aerial = list(
      comm_jaccard = pa_aerial,
      comm_bray    = comm_aerial,
      metadata     = metadata_aerial
    ),
    Aquatic = list(
      comm_jaccard = pa_aquatic,
      comm_bray    = comm_aquatic,
      metadata     = metadata_aquatic
    )
  )
  
  results <- list()
  counter <- 1
  total_tests <- 18
  
  for (comm_name in names(community_data)) {
    for (index_name in c("Jaccard", "Bray")) {
      for (component_name in c("Total", "Turnover", "Nestedness")) {
        
        cat("\n============================================================\n")
        cat("TEST", counter, "OF", total_tests, "\n")
        cat("Community:", comm_name,
            "| Index:", index_name,
            "| Component:", component_name, "\n")
        cat("============================================================\n")
        
        if (index_name == "Jaccard") {
          comm_matrix <- community_data[[comm_name]]$comm_jaccard
          method <- "jaccard"
        } else {
          comm_matrix <- community_data[[comm_name]]$comm_bray
          method <- "bray"
        }
        
        result <- permutation_spatial_beta(
          comm_matrix = comm_matrix,
          metadata = community_data[[comm_name]]$metadata,
          method = method,
          component = component_name,
          nperm = nperm,
          seed = seed + counter
        )
        
        results[[paste(comm_name, index_name, component_name,
                       sep = "_")]] <- result
        
        counter <- counter + 1
      }
    }
  }
  
  results
}

permutation_results <- run_all_spatial_permutations(
  nperm = 4999,
  seed = 123
)


################################################################################
## 7. SUMMARY TABLE  (raw p-values, no multiple-test correction)
################################################################################

permutation_summary <- bind_rows(
  lapply(
    names(permutation_results),
    function(x) {
      
      result <- permutation_results[[x]]
      parts <- strsplit(x, "_")[[1]]
      
      data.frame(
        Community     = parts[1],
        Index         = parts[2],
        Component     = paste(parts[-c(1, 2)], collapse = "_"),
        N_ERN         = result$n_ERN,
        N_IRN         = result$n_IRN,
        Mean_ERN      = result$observed_ERN,
        Mean_IRN      = result$observed_IRN,
        Difference_ERN_minus_IRN = result$observed_difference,
        P_two_sided   = result$p_value_two_sided,
        P_ERN_greater = result$p_value_ERN_greater,
        stringsAsFactors = FALSE
      )
    }
  )
)

permutation_summary_display <- permutation_summary %>%
  mutate(across(
    c(Mean_ERN, Mean_IRN, Difference_ERN_minus_IRN),
    ~ round(.x, 4)
  )) %>%
  mutate(across(
    c(P_two_sided, P_ERN_greater),
    ~ round(.x, 2)
  ))

cat("\n\n============================================================\n")
cat("SITE-LEVEL PERMUTATION TEST RESULTS (raw p-values)\n")
cat("============================================================\n\n")

print(permutation_summary_display)

write.csv(
  permutation_summary,
  "Spatial_Beta_Site_Level_Permutation_Results.csv",
  row.names = FALSE
)


################################################################################
## 8. SPATIAL BETA BY SAMPLING OCCASION
################################################################################

observed_period_results <- bind_rows(
  lapply(
    names(permutation_results),
    function(x) {
      
      parts <- strsplit(x, "_")[[1]]
      
      permutation_results[[x]]$observed_period %>%
        mutate(
          Community = parts[1],
          Index     = parts[2],
          Component = paste(parts[-c(1, 2)], collapse = "_")
        )
    }
  )
)

period_levels <- c("First time", "Second time", "Third time",
                   "Fourth time", "Fifth time", "Sixth time")

component_order <- c("Total", "Turnover", "Nestedness")
index_order     <- c("Bray", "Jaccard")
community_order <- c("Total", "Aerial", "Aquatic")
basin_order     <- c("ERN", "IRN")

observed_period_results <- observed_period_results %>%
  mutate(
    Period = factor(Period,
                    levels = period_levels,
                    labels = c("1st", "2nd", "3rd",
                               "4th", "5th", "6th")),
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order),
    Basin     = factor(Basin,     levels = basin_order)
  )

permutation_summary <- permutation_summary %>%
  mutate(
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order)
  )


################################################################################
## 9. TEMPORAL PLOT (diagnostic)
################################################################################

p_temporal <- ggplot(
  observed_period_results,
  aes(x = Period, y = Beta, group = Basin, color = Basin)
) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(
    aes(ymin = Beta - SE, ymax = Beta + SE),
    width = 0.12,
    linewidth = 0.6
  ) +
  geom_point(size = 2.5) +
  facet_grid(Component ~ Index + Community, scales = "free_y") +
  scale_color_manual(
    values = c("ERN" = "#D55E00", "IRN" = "#0072B2"),
    drop = FALSE
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Sampling occasion",
    y = expression("Mean spatial " * beta * " diversity (+- SE)"),
    color = "River basin"
  ) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold"),
    axis.text          = element_text(color = "black"),
    legend.position    = "bottom"
  )

print(p_temporal)


################################################################################
## 10. PUBLICATION BOXPLOT
##
## Descriptive only: inference comes from the site-level permutation
## test, not from treating pairwise values as independent.
################################################################################

spatial_p_labels <- permutation_summary %>%
  mutate(
    P_label = case_when(
      P_two_sided < 0.01 ~
        "p < 0.01",
      TRUE ~
        paste0("p = ", formatC(round(P_two_sided, 2),
                               format = "f", digits = 2))
    ),
    Significance = case_when(
      P_two_sided < 0.001 ~ "***",
      P_two_sided < 0.01  ~ "**",
      P_two_sided < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

## same vertical position for all p-values within a Component row
spatial_component_max <- observed_period_results %>%
  group_by(Component) %>%
  summarise(
    ymax = max(Beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    text_y         = ymax + 0.16 * abs(ymax),
    significance_y = ymax + 0.24 * abs(ymax)
  )

spatial_p_box_annot <- spatial_p_labels %>%
  select(Community, Index, Component, P_label, Significance) %>%
  distinct(Community, Index, Component, .keep_all = TRUE) %>%
  left_join(spatial_component_max, by = "Component") %>%
  mutate(
    x_text = 1.5,
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order)
  )

## panel letters (facet layout order):
##                BRAY                    JACCARD
##          Total Aerial Aquatic     Total Aerial Aquatic
## Total       A     B      C           D      E      F
## Turnover    G     H      I           J      K      L
## Nestedness  M     N      O           P      Q      R

spatial_panel_letters <- data.frame(
  Component = rep(component_order, each = 6),
  Index     = rep(rep(index_order, each = 3), times = 3),
  Community = rep(community_order, times = 6),
  Panel     = LETTERS[1:18],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order)
  )


################################################################################
## PLOT
################################################################################

p_spatial_beta <- ggplot(
  observed_period_results,
  aes(x = Basin, y = Beta, fill = Basin)
) +
  
  geom_boxplot(
    width = 0.55,
    alpha = 0.85,
    color = "black",
    linewidth = 0.6,
    outlier.shape = 21,
    outlier.size = 2.5,
    outlier.stroke = 0.5,
    outlier.fill = "white",
    outlier.alpha = 0.7
  ) +
  
  geom_jitter(width = 0.08, size = 2, alpha = 0.55) +
  
  geom_text(
    data = spatial_p_box_annot,
    aes(x = x_text, y = text_y, label = P_label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    color = "grey20"
  ) +
  
  geom_text(
    data = spatial_p_box_annot %>% filter(Significance != ""),
    aes(x = x_text, y = significance_y, label = Significance),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text(
    data = spatial_panel_letters,
    aes(x = -Inf, y = Inf, label = Panel),
    inherit.aes = FALSE,
    hjust = -0.25,
    vjust = 1.25,
    size = 5,
    fontface = "bold",
    color = "grey25"
  ) +
  
  facet_grid(
    Component ~ Index + Community,
    scales = "free_y",
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = c("ERN" = "#D55E00", "IRN" = "#0072B2"),
    drop = FALSE
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  
  labs(
    x = NULL,
    y = expression("Spatial " * beta * " diversity"),
    fill = "River basin"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA,
                                linewidth = 0.6),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90",
                                      linewidth = 0.3,
                                      linetype = "dashed"),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "black",
                                    linewidth = 0.6),
    axis.text  = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    panel.spacing   = unit(0.8, "lines"),
    plot.margin     = margin(15, 15, 15, 15)
  )

print(p_spatial_beta)


################################################################################
## 11. SAVE FIGURE
################################################################################

ggsave(
  filename = "FIGURE_4_Spatial_Beta_Basin_Permutation.tiff",
  plot = p_spatial_beta,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

cat("\n\nAnalysis completed.\n")

################################################################################
## END
################################################################################
