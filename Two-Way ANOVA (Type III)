############################################################
## ALPHA DIVERSITY + DENSITY — TYPE III ANOVA
## IMPROVED MODEL-BASED FIGURE (v2)
############################################################
#
# Same analysis pipeline as before.
# The figure section has been completely redesigned:
#
#   * facet_grid() guarantees aligned rows/columns
#   * row labels on the left, group names on top
#   * basins encoded by colour + shape + linetype
#   * panel letters (A-L) for publication
#   * BH-adjusted pairwise stars (ERN vs IRN) drawn on panels
#   * single merged legend, light y-grid, styled strips
#
############################################################

rm(list = ls())

############################################################
## 0. PACKAGES
############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(car)
library(emmeans)
library(openxlsx)
library(grid)


############################################################
## 1. READ DATA
############################################################

file <- "FUNCAP_FWB.xlsx"

df_raw <- read_excel(
  file,
  sheet = "macroinvertebrates",
  col_names = FALSE
)


############################################################
## 2. DEFINE COLUMN NAMES
############################################################

colnames(df_raw) <- as.character(df_raw[4, ])
df <- df_raw[-c(1:4), ]


############################################################
## 3. METADATA AND TAXA COLUMNS
############################################################

meta_cols <- c("SU", "Basin", "Flow", "Period", "Date", "Codes")

taxa_cols <- colnames(df)[!colnames(df) %in% meta_cols]


############################################################
## 4. FUNCTIONAL GROUPS
############################################################

fg_df <- data.frame(
  taxa  = taxa_cols,
  group = as.character(
    df_raw[3, match(taxa_cols, colnames(df_raw))]
  ),
  stringsAsFactors = FALSE
)

aerial_taxa  <- fg_df$taxa[fg_df$group == "Aerial"]
aquatic_taxa <- fg_df$taxa[fg_df$group == "Purely Aquatic"]


############################################################
## 5. TAXA TO NUMERIC
############################################################

df[taxa_cols] <- lapply(df[taxa_cols], as.numeric)


############################################################
## 6. GROUP DATASETS
############################################################

macro_total   <- df[, c(meta_cols, taxa_cols)]
macro_aerial  <- df[, c(meta_cols, aerial_taxa)]
macro_aquatic <- df[, c(meta_cols, aquatic_taxa)]

metadata     <- macro_total[, meta_cols]
comm_total   <- macro_total[, taxa_cols]
comm_aerial  <- macro_aerial[, aerial_taxa]
comm_aquatic <- macro_aquatic[, aquatic_taxa]


############################################################
## 7. REMOVE EMPTY SAMPLING UNITS
############################################################

all_na   <- rowSums(is.na(comm_total)) == ncol(comm_total)
all_zero <- rowSums(comm_total, na.rm = TRUE) == 0
keep     <- !(all_na | all_zero)

metadata     <- metadata[keep, , drop = FALSE]
macro_total  <- macro_total[keep, , drop = FALSE]
macro_aerial <- macro_aerial[keep, , drop = FALSE]
macro_aquatic <- macro_aquatic[keep, , drop = FALSE]
comm_total   <- comm_total[keep, , drop = FALSE]
comm_aerial  <- comm_aerial[keep, , drop = FALSE]
comm_aquatic <- comm_aquatic[keep, , drop = FALSE]


############################################################
## 8. METRICS
############################################################

calc_metrics <- function(comm_matrix, metadata = NULL) {
  
  Density <- rowSums(comm_matrix, na.rm = TRUE)
  q0      <- specnumber(comm_matrix)
  H       <- diversity(comm_matrix, index = "shannon")
  q1      <- exp(H)
  D       <- diversity(comm_matrix, index = "simpson")
  q2      <- 1 / (1 - D)
  
  metrics <- data.frame(Density = Density, q0 = q0, q1 = q1, q2 = q2)
  
  if (!is.null(metadata)) metrics <- cbind(metadata, metrics)
  metrics
}

metrics_total   <- calc_metrics(comm_total,   metadata)
metrics_aerial  <- calc_metrics(comm_aerial,  metadata)
metrics_aquatic <- calc_metrics(comm_aquatic, metadata)

metrics_total$Group   <- "Total"
metrics_aerial$Group  <- "Aerial"
metrics_aquatic$Group <- "Purely Aquatic"

metrics_all <- bind_rows(metrics_total, metrics_aerial, metrics_aquatic)


############################################################
## 9. COMMON PERIODS + ANALYSIS DATASETS
############################################################

common_periods <- c("First time", "Fifth time", "Sixth time")

metrics_common <- metrics_all %>%
  filter(Period %in% common_periods) %>%
  mutate(
    Basin  = factor(Basin,  levels = c("ERN", "IRN")),
    Group  = factor(Group,  levels = c("Aerial", "Purely Aquatic", "Total")),
    Period = factor(Period, levels = common_periods)
  )

aquatic_exclude_su <- c("BN07-C5", "BN10-C5")

total_dat <- metrics_common %>%
  filter(Group == "Total") %>% droplevels()

aerial_dat <- metrics_common %>%
  filter(Group == "Aerial") %>% droplevels()

aquatic_dat <- metrics_common %>%
  filter(Group == "Purely Aquatic") %>%
  filter(!SU %in% aquatic_exclude_su) %>%
  droplevels()

stopifnot(
  nrow(aquatic_dat %>% filter(SU %in% aquatic_exclude_su)) == 0
)


############################################################
## 10. TYPE III CONTRASTS
############################################################

options(contrasts = c("contr.sum", "contr.poly"))


############################################################
## 11. FIT MODELS
############################################################

fit_group <- function(dat, group_name) {
  
  models <- list(
    Density = lm(log1p(Density) ~ Basin * Period, data = dat),
    q0      = lm(q0  ~ Basin * Period, data = dat),
    q1      = lm(q1  ~ Basin * Period, data = dat),
    q2      = lm(q2  ~ Basin * Period, data = dat)
  )
  
  anovas <- lapply(models, function(m) car::Anova(m, type = "III"))
  
  list(group = group_name, models = models, anova = anovas)
}

results_total   <- fit_group(total_dat,   "Total")
results_aerial  <- fit_group(aerial_dat,  "Aerial")
results_aquatic <- fit_group(aquatic_dat, "Purely Aquatic")


############################################################
## 12. PRINT ANOVA TABLES
############################################################

for (res in list(results_total, results_aerial, results_aquatic)) {
  
  cat("\n\n========================================================\n")
  cat("TYPE III ANOVA —", res$group, "\n")
  cat("========================================================\n")
  
  for (metric in c("Density", "q0", "q1", "q2")) {
    cat("\n---", metric, "---\n")
    print(res$anova[[metric]])
  }
}


############################################################
## 13. ASSUMPTIONS
############################################################

check_assumptions <- function(model, dat, response_name, group_name) {
  
  res      <- residuals(model)
  shapiro  <- shapiro.test(res)
  
  diagnostic_dat <- dat %>%
    mutate(
      Residuals = res,
      Grouping  = interaction(Basin, Period, drop = TRUE)
    )
  
  bf <- car::leveneTest(Residuals ~ Grouping, data = diagnostic_dat, center = median)
  
  cook      <- cooks.distance(model)
  leverage  <- hatvalues(model)
  std_res   <- rstandard(model)
  
  data.frame(
    Group              = group_name,
    Metric             = response_name,
    N                  = nrow(dat),
    Shapiro_W          = unname(shapiro$statistic),
    Shapiro_P          = shapiro$p.value,
    Brown_Forsythe_F   = bf$`F value`[1],
    Brown_Forsythe_P   = bf$`Pr(>F)`[1],
    High_Cook_N        = sum(cook > 4 / nrow(dat), na.rm = TRUE),
    High_Leverage_N    = sum(leverage > (2 * length(coef(model))) / nrow(dat), na.rm = TRUE),
    High_Residual_N    = sum(abs(std_res) > 3, na.rm = TRUE)
  )
}

assumptions <- bind_rows(
  check_assumptions(results_total$models$Density,   total_dat,   "Density", "Total"),
  check_assumptions(results_total$models$q0,        total_dat,   "q0",      "Total"),
  check_assumptions(results_total$models$q1,        total_dat,   "q1",      "Total"),
  check_assumptions(results_total$models$q2,        total_dat,   "q2",      "Total"),
  check_assumptions(results_aerial$models$Density,  aerial_dat,  "Density", "Aerial"),
  check_assumptions(results_aerial$models$q0,       aerial_dat,  "q0",      "Aerial"),
  check_assumptions(results_aerial$models$q1,       aerial_dat,  "q1",      "Aerial"),
  check_assumptions(results_aerial$models$q2,       aerial_dat,  "q2",      "Aerial"),
  check_assumptions(results_aquatic$models$Density, aquatic_dat, "Density", "Purely Aquatic"),
  check_assumptions(results_aquatic$models$q0,      aquatic_dat, "q0",      "Purely Aquatic"),
  check_assumptions(results_aquatic$models$q1,      aquatic_dat, "q1",      "Purely Aquatic"),
  check_assumptions(results_aquatic$models$q2,      aquatic_dat, "q2",      "Purely Aquatic")
)

print(assumptions)


############################################################
## 14. ESTIMATED MARGINAL MEANS (BACK-TRANSFORMED)
############################################################

period_short <- c(
  "First time" = "First",
  "Fifth time" = "Fifth",
  "Sixth time" = "Sixth"
)

get_emm <- function(model, metric, group) {
  
  emmeans(model, ~ Basin * Period) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    mutate(Metric = metric, Group = group)
}

emm_all <- bind_rows(
  get_emm(results_total$models$Density,   "Density", "Total"),
  get_emm(results_total$models$q0,        "q0",      "Total"),
  get_emm(results_total$models$q1,        "q1",      "Total"),
  get_emm(results_total$models$q2,        "q2",      "Total"),
  get_emm(results_aerial$models$Density,  "Density", "Aerial"),
  get_emm(results_aerial$models$q0,       "q0",      "Aerial"),
  get_emm(results_aerial$models$q1,       "q1",      "Aerial"),
  get_emm(results_aerial$models$q2,       "q2",      "Aerial"),
  get_emm(results_aquatic$models$Density, "Density", "Purely Aquatic"),
  get_emm(results_aquatic$models$q0,      "q0",      "Purely Aquatic"),
  get_emm(results_aquatic$models$q1,      "q1",      "Purely Aquatic"),
  get_emm(results_aquatic$models$q2,      "q2",      "Purely Aquatic")
) %>%
  mutate(
    response       = if_else(Metric == "Density", expm1(emmean),   emmean),
    lower_response = if_else(Metric == "Density", expm1(lower.CL), lower.CL),
    upper_response = if_else(Metric == "Density", expm1(upper.CL), upper.CL),
    Group  = factor(Group,  levels = c("Aerial", "Purely Aquatic", "Total")),
    Metric = factor(Metric, levels = c("Density", "q0", "q1", "q2")),
    Period = factor(period_short[as.character(Period)],
                    levels = unname(period_short)),
    Basin  = factor(as.character(Basin), levels = c("ERN", "IRN"))
  )

stopifnot(nrow(emm_all) == 72)


############################################################
## 15. PAIRWISE COMPARISONS (ERN vs IRN WITHIN PERIOD)
############################################################

get_pw <- function(model, metric, group) {
  
  emmeans(model, ~ Basin | Period) %>%
    pairs(adjust = "BH") %>%
    summary() %>%
    as.data.frame() %>%
    mutate(Metric = metric, Group = group)
}

pairwise_all <- bind_rows(
  get_pw(results_total$models$Density,   "Density", "Total"),
  get_pw(results_total$models$q0,        "q0",      "Total"),
  get_pw(results_total$models$q1,        "q1",      "Total"),
  get_pw(results_total$models$q2,        "q2",      "Total"),
  get_pw(results_aerial$models$Density,  "Density", "Aerial"),
  get_pw(results_aerial$models$q0,       "q0",      "Aerial"),
  get_pw(results_aerial$models$q1,       "q1",      "Aerial"),
  get_pw(results_aerial$models$q2,       "q2",      "Aerial"),
  get_pw(results_aquatic$models$Density, "Density", "Purely Aquatic"),
  get_pw(results_aquatic$models$q0,      "q0",      "Purely Aquatic"),
  get_pw(results_aquatic$models$q1,      "q1",      "Purely Aquatic"),
  get_pw(results_aquatic$models$q2,      "q2",      "Purely Aquatic")
) %>%
  mutate(
    Group  = factor(Group,  levels = c("Aerial", "Purely Aquatic", "Total")),
    Metric = factor(Metric, levels = c("Density", "q0", "q1", "q2")),
    Period = factor(period_short[as.character(Period)],
                    levels = unname(period_short)),
    Significance = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      p.value <= 0.10  ~ ".",
      TRUE ~ ""
    )
  )

print(pairwise_all)


############################################################
## 16. SIGNIFICANCE STAR POSITIONS
############################################################

sig_dat <- pairwise_all %>%
  filter(Significance != "") %>%
  left_join(
    emm_all %>%
      group_by(Group, Metric, Period) %>%
      summarise(ymax = max(upper_response), .groups = "drop"),
    by = c("Group", "Metric", "Period")
  ) %>%
  mutate(y = ymax * 1.05)


############################################################
## 16.1 TYPE III ANOVA EFFECTS FOR PANEL ANNOTATION
############################################################
#
# Extracts the p-values of the three focal Type III effects
# (Basin, Period, Basin x Period) for every Group x Metric
# panel and formats them as a compact annotation label.
#
############################################################

format_stars <- function(p) {
  case_when(
    is.na(p)    ~ "",
    p <= 0.001  ~ "***",
    p <= 0.01   ~ "**",
    p <= 0.05   ~ "*",
    p <= 0.10   ~ ".",
    TRUE        ~ "ns"
  )
}

extract_effects <- function(res) {
  
  bind_rows(
    lapply(
      names(res$anova),
      function(metric) {
        
        a <- as.data.frame(
          res$anova[[metric]]
        )
        
        a$Effect <- rownames(a)
        rownames(a) <- NULL
        
        a %>%
          transmute(
            Metric  = metric,
            Effect,
            p_value = `Pr(>F)`
          )
      }
    )
  )
}

anova_effects <- bind_rows(
  
  extract_effects(results_total) %>%
    mutate(Group = "Total"),
  
  extract_effects(results_aerial) %>%
    mutate(Group = "Aerial"),
  
  extract_effects(results_aquatic) %>%
    mutate(Group = "Purely Aquatic")
  
) %>%
  
  filter(
    Effect %in% c(
      "Basin",
      "Period",
      "Basin:Period"
    )
  ) %>%
  
  mutate(
    
    Group = factor(
      Group,
      levels = c(
        "Aerial",
        "Purely Aquatic",
        "Total"
      )
    ),
    
    Metric = factor(
      Metric,
      levels = c(
        "Density",
        "q0",
        "q1",
        "q2"
      )
    ),
    
    term = case_when(
      Effect == "Basin"        ~ "Basin",
      Effect == "Period"       ~ "Period",
      Effect == "Basin:Period" ~ "Basin x Period"
    ),
    
    sig = format_stars(p_value)
    
  ) %>%
  
  ## one short label per term -------------------------------
mutate(
  label = paste0(term, " ", sig)
) %>%
  
  ## flag panels with at least one significant effect -------
group_by(
  Group,
  Metric
) %>%
  
  mutate(
    panel_significant = any(
      sig %in% c("*", "**", "***", ".")
    )
  ) %>%
  
  ## collapse all three terms into one label per panel ------
summarise(
  label = paste(
    label,
    collapse = "\n"
  ),
  panel_significant = first(
    panel_significant
  ),
  .groups = "drop"
) %>%
  
  ## keep only panels with a significant effect -------------
## (C, H and K with the current results)
filter(
  panel_significant
) %>%
  
  select(
    -panel_significant
  )

print(anova_effects)


############################################################
## 17. PANEL LETTERS (A-L)
############################################################

panel_letters <- data.frame(
  Metric = rep(levels(emm_all$Metric), each = 3),
  Group  = rep(levels(emm_all$Group),  times = 4),
  lab    = LETTERS[1:12],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Metric = factor(Metric, levels = levels(emm_all$Metric)),
    Group  = factor(Group,  levels = levels(emm_all$Group))
  )


############################################################
## 18. ROW LABELS (PARSED EXPRESSIONS)
############################################################

metric_labeller <- as_labeller(
  c(
    Density = "Density~(ind.~m^{-2})",
    q0      = "Species~richness~(q[0])",
    q1      = "Effective~Shannon~diversity~(q[1])",
    q2      = "Inverse~Simpson~diversity~(q[2])"
  ),
  default = label_parsed
)


############################################################
## 19. FIGURE 3 - FRESHWATER BIOLOGY
############################################################

basin_pal <- c(ERN = "#D55E00", IRN = "#0072B2")

figure_alpha_v2 <- ggplot(
  emm_all,
  aes(
    x = Period,
    y = response,
    colour = Basin,
    shape  = Basin,
    linetype = Basin,
    group  = Basin
  )
) +
  
  ## 95% CI ------------------------------------------------
geom_errorbar(
  aes(ymin = lower_response, ymax = upper_response),
  width = 0.09,
  linewidth = 0.55,
  alpha = 0.9,
  show.legend = FALSE
) +
  
  ## trajectories -------------------------------------------
geom_line(linewidth = 0.7) +
  
  ## model estimates ----------------------------------------
geom_point(size = 2.8, stroke = 0.75) +
  
  ## significance stars -------------------------------------
geom_text(
  data = sig_dat,
  aes(x = Period, y = y, label = Significance),
  inherit.aes = FALSE,
  size = 3.6,
  fontface = "bold",
  colour = "grey15",
  vjust = 0
) +
  
  ## panel letters ------------------------------------------
geom_text(
  data = panel_letters,
  aes(x = -Inf, y = Inf, label = lab),
  inherit.aes = FALSE,
  hjust = -0.5,
  vjust = 1.4,
  fontface = "bold",
  size = 3.8,
  colour = "black"
) +
  
  ## Type III two-way ANOVA effects -------------------------
geom_text(
  data = anova_effects,
  aes(
    x = Inf,
    y = Inf,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = 1.15,
  vjust = 1.5,
  size = 2.7,
  lineheight = 0.95,
  colour = "grey25"
) +
  
  ## scales -------------------------------------------------
scale_colour_manual(
  values = basin_pal,
  name = "River basin"
) +
  scale_shape_manual(
    values = c(ERN = 21, IRN = 22),
    name = "River basin"
  ) +
  scale_linetype_manual(
    values = c(ERN = "solid", IRN = "longdash"),
    name = "River basin"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.30))) +
  
  ## layout: rows = metrics, columns = groups ---------------
facet_grid(
  Metric ~ Group,
  scales   = "free_y",
  switch   = "y",
  labeller = labeller(Metric = metric_labeller)
) +
  
  ## labels -------------------------------------------------
labs(
  x = "Sampling period",
  y = NULL#,
  #caption = paste(
   # "Points are estimated marginal means; bars are 95% CI.",
    #"Stars above periods: BH-adjusted ERN vs IRN pairwise tests.",
    #"Text at top right of panels C, H and K: full Type III two-way ANOVA results",
    #"(Basin, Period and Basin x Period;",
   # "* p < 0.05, ** p < 0.01, *** p < 0.001, ns = not significant)."
#  )
) +
  
  ## theme --------------------------------------------------
theme_classic(base_size = 11) +
  theme(
    axis.title.x      = element_text(face = "bold", size = 11,
                                     margin = margin(t = 6)),
    axis.text         = element_text(colour = "grey20", size = 9.5),
    axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line         = element_line(linewidth = 0.4, colour = "grey40"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.4),
    panel.spacing.x   = unit(1.0, "lines"),
    panel.spacing.y   = unit(0.8, "lines"),
    strip.background  = element_rect(fill = "grey95", colour = "grey65",
                                     linewidth = 0.35),
    strip.text        = element_text(face = "bold", size = 10.5,
                                     colour = "grey15"),
    strip.text.y.left = element_text(face = "bold", size = 10.5,
                                     colour = "grey15", angle = 0, hjust = 1),
    strip.placement   = "outside",
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9.5),
    legend.key.width  = unit(1.3, "cm"),
    legend.key.height = unit(0.45, "cm"),
    plot.caption      = element_text(size = 8, colour = "grey35",
                                     hjust = 0, margin = margin(t = 8)),
    plot.margin       = margin(6, 10, 6, 10)
  )

print(figure_alpha_v2)


############################################################
## 20. SAVE FIGURE
############################################################
ggsave(
  filename = "Alpha_Diversity_Figure_v2.tiff",
  plot     = figure_alpha_v2,
  width    = 9.3,
  height   = 9.8,
  units    = "in",
  dpi      = 600,
  compression = "lzw"
)


############################################################
## 21. EXPORT STATISTICS (SAME AS BEFORE)
############################################################

extract_anova <- function(anova_object, group_name, metric_name) {
  
  out <- as.data.frame(anova_object)
  out$Effect <- rownames(out)
  rownames(out) <- NULL
  out$Group  <- group_name
  out$Metric <- metric_name
  out
}

anova_all <- bind_rows(
  extract_anova(results_total$anova$Density,   "Total",            "Density"),
  extract_anova(results_total$anova$q0,        "Total",            "q0"),
  extract_anova(results_total$anova$q1,        "Total",            "q1"),
  extract_anova(results_total$anova$q2,        "Total",            "q2"),
  extract_anova(results_aerial$anova$Density,  "Aerial",           "Density"),
  extract_anova(results_aerial$anova$q0,       "Aerial",           "q0"),
  extract_anova(results_aerial$anova$q1,       "Aerial",           "q1"),
  extract_anova(results_aerial$anova$q2,       "Aerial",           "q2"),
  extract_anova(results_aquatic$anova$Density, "Purely Aquatic",   "Density"),
  extract_anova(results_aquatic$anova$q0,      "Purely Aquatic",   "q0"),
  extract_anova(results_aquatic$anova$q1,      "Purely Aquatic",   "q1"),
  extract_anova(results_aquatic$anova$q2,      "Purely Aquatic",   "q2")
)

write.xlsx(
  list(
    Type_III_ANOVA             = anova_all,
    Assumptions                = assumptions,
    Pairwise_EMMeans           = pairwise_all,
    Model_EMMeans              = emm_all,
    Analysis_Data_Total        = total_dat,
    Analysis_Data_Aerial       = aerial_dat,
    Analysis_Data_Purely_Aquatic = aquatic_dat
  ),
  file      = "Alpha_Diversity_TypeIII_ANOVA_v2.xlsx",
  overwrite = TRUE
)

############################################################
## END
############################################################
