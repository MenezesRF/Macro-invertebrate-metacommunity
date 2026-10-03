## MENEZES et al. (Freshwater Biology manuscript)
# The scripts below should be run after the scripts in the file "P2_spatial_beta_permutation_test.R".

################################################################################
##
##  DESCRIPTION:
##  Codes for testing P2 and generating FIGURE 5, that spatial β-diversity is higher ##  in ERN than in IRN whereas temporal β-diversity is lower in ERN than in IRN
##  
################################################################################
##
##  TEMPORAL BETA DIVERSITY — ERN vs IRN
##
##  Temporal beta is calculated for the SAME SITE across CONSECUTIVE
##  sampling occasions (1st->2nd, 2nd->3rd, ...). Sites do NOT need
##  all six occasions.
##
##  Statistical unit = the SITE: each site contributes ONE value
##  (mean of its transitions) to the basin comparison.
##
##  Inference: site-level permutation test (4,999 permutations),
##  randomizing basin labels among sites while keeping n_ERN = 9 and
##  n_IRN = 18. reported p-values are raw two-sided permutation p-values.
##
##  REQUIRES objects from the data-preparation step:
##    comm_total / comm_aerial / comm_aquatic   (abundance matrices)
##    pa_total / pa_aerial / pa_aquatic         (presence-absence)
##    metadata_total / metadata_aerial / metadata_aquatic
##  (e.g., run the simplified spatial script first, WITHOUT its
##   rm(list = ls()) line.)
## Due to the time it takes to run the analyses, it is recommended that a lower number of permutations (nperm) be used for the replication of the analyses. 
################################################################################


################################################################################
## 1. TEMPORAL BETA FOR ONE SITE'S CONSECUTIVE OCCASION PAIRS
################################################################################

calculate_temporal_beta <- function(
    comm_matrix,
    metadata,
    method = "jaccard",
    component = "Total"
) {
  
  if (!"Site_ID" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(Site_ID = sub("-C[0-9]+$", "", as.character(SU)))
  }
  
  period_order <- c(
    "First time", "Second time", "Third time",
    "Fourth time", "Fifth time", "Sixth time"
  )
  
  data_work <- metadata %>%
    mutate(
      Row_ID = seq_len(n()),
      Period = factor(Period, levels = period_order)
    )
  
  results <- list()
  
  for (site in unique(data_work$Site_ID)) {
    
    site_data <- data_work %>%
      filter(Site_ID == site) %>%
      arrange(Period)
    
    if (nrow(site_data) < 2) next
    
    for (i in 1:(nrow(site_data) - 1)) {
      
      row1 <- site_data$Row_ID[i]
      row2 <- site_data$Row_ID[i + 1]
      
      period1 <- as.character(site_data$Period[i])
      period2 <- as.character(site_data$Period[i + 1])
      
      comm_pair <- comm_matrix[c(row1, row2), , drop = FALSE]
      
      keep_species <- colSums(comm_pair, na.rm = TRUE) > 0
      comm_pair <- comm_pair[, keep_species, drop = FALSE]
      
      if (ncol(comm_pair) < 2) next
      
      if (method == "jaccard") {
        
        comm_pair <- decostand(comm_pair, method = "pa")
        
        beta <- beta.pair(comm_pair, index.family = "jaccard")
        
        beta_value <- switch(
          component,
          Total      = as.numeric(beta$beta.jac[1]),
          Turnover   = as.numeric(beta$beta.jtu[1]),
          Nestedness = as.numeric(beta$beta.jne[1])
        )
      }
      
      if (method == "bray") {
        
        beta <- beta.pair.abund(comm_pair, index.family = "bray")
        
        beta_value <- switch(
          component,
          Total      = as.numeric(beta$beta.bray[1]),
          Turnover   = as.numeric(beta$beta.bray.bal[1]),
          Nestedness = as.numeric(beta$beta.bray.gra[1])
        )
      }
      
      results[[length(results) + 1]] <- data.frame(
        Site_ID    = site,
        Basin      = as.character(site_data$Basin[i]),
        Period1    = period1,
        Period2    = period2,
        Transition = paste0(period1, " -> ", period2),
        Beta       = beta_value,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(results)
}


################################################################################
## 2. ALL THREE COMPONENTS FOR ONE COMMUNITY x INDEX COMBINATION
################################################################################

calculate_all_temporal <- function(
    comm_matrix,
    metadata,
    community_name,
    method,
    pa_matrix = NULL
) {
  
  components <- c("Total", "Turnover", "Nestedness")
  results <- list()
  
  for (comp in components) {
    
    if (method == "jaccard") {
      result <- calculate_temporal_beta(
        pa_matrix, metadata,
        method = "jaccard", component = comp
      )
    } else {
      result <- calculate_temporal_beta(
        comm_matrix, metadata,
        method = "bray", component = comp
      )
    }
    
    result$Community <- community_name
    result$Index     <- ifelse(method == "jaccard", "Jaccard", "Bray")
    result$Component <- comp
    
    results[[comp]] <- result
  }
  
  bind_rows(results)
}


################################################################################
## 3. ALL TEMPORAL BETA RESULTS (6 community x index runs)
################################################################################

temporal_all <- bind_rows(
  
  calculate_all_temporal(comm_total,   metadata_total,
                         "Total",   "jaccard", pa_total),
  calculate_all_temporal(comm_total,   metadata_total,
                         "Total",   "bray"),
  calculate_all_temporal(comm_aerial,  metadata_aerial,
                         "Aerial",  "jaccard", pa_aerial),
  calculate_all_temporal(comm_aerial,  metadata_aerial,
                         "Aerial",  "bray"),
  calculate_all_temporal(comm_aquatic, metadata_aquatic,
                         "Aquatic", "jaccard", pa_aquatic),
  calculate_all_temporal(comm_aquatic, metadata_aquatic,
                         "Aquatic", "bray")
)


################################################################################
## 4. ONE VALUE PER SITE
################################################################################

temporal_site <- temporal_all %>%
  group_by(Site_ID, Basin, Community, Index, Component) %>%
  summarise(
    Mean_Temporal_Beta = mean(Beta, na.rm = TRUE),
    SD_Temporal_Beta   = sd(Beta, na.rm = TRUE),
    N_Transitions      = sum(!is.na(Beta)),
    .groups = "drop"
  )


################################################################################
## 5. SITE-LEVEL PERMUTATION TEST
################################################################################

permutation_temporal_beta <- function(data, nperm = 4999, seed = 123) {
  
  set.seed(seed)
  
  observed <- data %>%
    group_by(Basin) %>%
    summarise(
      Mean = mean(Mean_Temporal_Beta, na.rm = TRUE),
      .groups = "drop"
    )
  
  observed_ERN <- observed$Mean[observed$Basin == "ERN"]
  observed_IRN <- observed$Mean[observed$Basin == "IRN"]
  
  observed_difference <- observed_ERN - observed_IRN
  
  sites <- data %>% distinct(Site_ID, Basin)
  
  n_ERN <- sum(sites$Basin == "ERN")
  n_IRN <- sum(sites$Basin == "IRN")
  
  null_difference <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    
    shuffled_sites <- sample(sites$Site_ID)
    
    pseudo_ERN <- shuffled_sites[1:n_ERN]
    
    perm_means <- data %>%
      mutate(
        Perm_Basin = ifelse(Site_ID %in% pseudo_ERN, "ERN", "IRN")
      ) %>%
      group_by(Perm_Basin) %>%
      summarise(
        Mean = mean(Mean_Temporal_Beta, na.rm = TRUE),
        .groups = "drop"
      )
    
    perm_ERN <- perm_means$Mean[perm_means$Perm_Basin == "ERN"]
    perm_IRN <- perm_means$Mean[perm_means$Perm_Basin == "IRN"]
    
    null_difference[i] <- perm_ERN - perm_IRN
  }
  
  p_two_sided <- (
    sum(abs(null_difference) >= abs(observed_difference)) + 1
  ) / (nperm + 1)
  
  list(
    observed_ERN        = observed_ERN,
    observed_IRN        = observed_IRN,
    observed_difference = observed_difference,
    null_difference     = null_difference,
    p_two_sided         = p_two_sided,
    n_ERN               = n_ERN,
    n_IRN               = n_IRN
  )
}


################################################################################
## 6. RUN ALL 18 TESTS (3 communities x 2 indices x 3 components)
################################################################################

temporal_combinations <- temporal_site %>%
  distinct(Community, Index, Component)

temporal_permutation_results <- list()

for (i in seq_len(nrow(temporal_combinations))) {
  
  comm <- temporal_combinations$Community[i]
  idx  <- temporal_combinations$Index[i]
  comp <- temporal_combinations$Component[i]
  
  data_subset <- temporal_site %>%
    filter(Community == comm, Index == idx, Component == comp)
  
  key <- paste(comm, idx, comp, sep = "_")
  
  cat("\nRunning permutation:", key, "\n")
  
  temporal_permutation_results[[key]] <-
    permutation_temporal_beta(
      data_subset,
      nperm = 4999,
      seed = 123
    )
}


################################################################################
## 7. SUMMARY TABLE
################################################################################

temporal_permutation_summary <- bind_rows(
  lapply(
    names(temporal_permutation_results),
    function(key) {
      
      res <- temporal_permutation_results[[key]]
      parts <- strsplit(key, "_")[[1]]
      
      data.frame(
        Community                = parts[1],
        Index                    = parts[2],
        Component                = parts[3],
        Mean_ERN                 = res$observed_ERN,
        Mean_IRN                 = res$observed_IRN,
        Difference_ERN_minus_IRN = res$observed_difference,
        P_value                  = res$p_two_sided,
        N_ERN                    = res$n_ERN,
        N_IRN                    = res$n_IRN,
        stringsAsFactors = FALSE
      )
    }
  )
)

print(temporal_permutation_summary)

write.csv(
  temporal_permutation_summary,
  "Temporal_Beta_Site_Level_Permutation_Results.csv",
  row.names = FALSE
)


################################################################################
## 8. FACTOR ORDERS
################################################################################

component_order <- c("Total", "Turnover", "Nestedness")
index_order     <- c("Bray", "Jaccard")
community_order <- c("Total", "Aerial", "Aquatic")
basin_order     <- c("ERN", "IRN")

temporal_beta_long <- temporal_site %>%
  rename(Beta = Mean_Temporal_Beta) %>%
  mutate(
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order),
    Basin     = factor(Basin,     levels = basin_order)
  )

temporal_permutation_summary <- temporal_permutation_summary %>%
  mutate(
    Component = factor(Component, levels = component_order),
    Index     = factor(Index,     levels = index_order),
    Community = factor(Community, levels = community_order)
  )


################################################################################
## 9. P-VALUE AND SIGNIFICANCE LABELS
################################################################################

temporal_p_labels <- temporal_permutation_summary %>%
  mutate(
    P_label = case_when(
      P_value < 0.001 ~
        "p < 0.001",
      TRUE ~
        paste0("p = ", formatC(P_value, format = "f", digits = 3))
    ),
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )


################################################################################
## 10. ANNOTATION POSITIONS AND PANEL LETTERS
################################################################################

## same vertical position for all p-values within a Component row
temporal_p_box_annot <- temporal_beta_long %>%
  group_by(Component) %>%
  summarise(
    max_beta = max(Beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_position     = max_beta + 0.15 * abs(max_beta),
    y_significance = max_beta + 0.23 * abs(max_beta)
  ) %>%
  right_join(
    temporal_p_labels %>%
      select(Community, Index, Component, P_label, Significance),
    by = "Component"
  ) %>%
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

temporal_panel_letters <- data.frame(
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
## 11. PUBLICATION BOXPLOT
################################################################################

p_temporal_box <- ggplot(
  temporal_beta_long,
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
  
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  
  geom_text(
    data = temporal_p_box_annot,
    aes(x = x_text, y = y_position, label = P_label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    color = "grey20"
  ) +
  
  geom_text(
    data = temporal_p_box_annot %>% filter(Significance != ""),
    aes(x = x_text, y = y_significance, label = Significance),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  geom_text(
    data = temporal_panel_letters,
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
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.38))
  ) +
  
  scale_fill_manual(
    values = c("ERN" = "#D55E00", "IRN" = "#0072B2"),
    drop = FALSE
  ) +
  
  labs(
    x = NULL,
    y = expression("Temporal " * beta * " diversity"),
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

print(p_temporal_box)


################################################################################
## 12. SAVE FIGURE
################################################################################
ggsave(
  filename = "FIGURE_5_Temporal_Beta_Basin_Permutation.tiff",
  plot = p_temporal_box,
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
