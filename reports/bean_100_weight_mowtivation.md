Soybean 100 seed weight
================

- [Setup](#setup)
- [Packages](#packages)
- [Data import and preparation](#data-import-and-preparation)
- [Model testing and diagnostics](#model-testing-and-diagnostics)
  - [Exploratory checks](#exploratory-checks)
  - [Model selection](#model-selection)
  - [Diagnostics (residuals, dispersion, zero
    inflation)](#diagnostics-residuals-dispersion-zero-inflation)
- [Summary tables](#summary-tables)
  - [Post-hoc treatment means and
    CLDs](#post-hoc-treatment-means-and-clds)
  - [Global response summary (ANOVA, AIC, model
    info)](#global-response-summary-anova-aic-model-info)
- [Figures](#figures)
  - [Model-predicted means (g 100 seeds
    )](#model-predicted-means-g-100-seeds-)
  - [Raw means (g 100 seeds )](#raw-means-g-100-seeds-)

# Setup

# Packages

``` r
# Packages
library(tidyverse)    # includes dplyr, ggplot2, readr, tibble, etc.
library(janitor)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(multcomp)
library(car)
library(kableExtra)
library(here)
library(conflicted)
library(lme4)
library(WrensBookshelf)
library(writexl)

# Handle conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

# Treatment level order (use everywhere)
mow_levels <- c(
  "Rolled, no control",
  "Rolled, mowing",
  "Rolled, high-residue cultivation",
  "Tilled, mowing",
  "Tilled, cultivation"
)

# One consistent CVD-safe color palette for all figures (WrensBookshelf)
fill_cols <- WB_brewer(
  name = "WhatWellBuild",
  n    = length(mow_levels),
  type = "discrete"
) |>
  setNames(mow_levels)

# x-axis label helpers ---------------------------------------------------

# Break on spaces (if you ever want every word on its own line)
label_break_spaces <- function(x) {
  stringr::str_replace_all(x, " ", "\n")
}

# Break after the comma: "Rolled,\nno control", etc.
label_break_comma <- function(x) {
  stringr::str_replace_all(x, ", ", ",\n")
}

# Break after comma AND split "high-residue cultivation"
# -> "Rolled,\nhigh-residue\ncultivation"
label_break_comma_cult <- function(x) {
  x |>
    stringr::str_replace("high-residue cultivation",
                         "high-residue\ncultivation") |>
    stringr::str_replace_all(", ", ",\n")
}

# Helper: tidy emmeans output regardless of CI column names --------------
# (works directly on an emmeans object)

tidy_emm <- function(emm, ref_levels = NULL) {
  emm_df <- as.data.frame(emm)

  lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(emm_df))[1]
  ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(emm_df))[1]

  if (is.na(lcl_col) || is.na(ucl_col)) {
    stop("Could not find CI columns in emmeans output.")
  }

  out <- emm_df |>
    dplyr::mutate(
      ci_low  = .data[[lcl_col]],
      ci_high = .data[[ucl_col]]
    )

  if (!is.null(ref_levels) && "weed_trt" %in% names(out)) {
    out <- out |>
      dplyr::mutate(weed_trt = factor(weed_trt, levels = ref_levels))
  }

  out
}
```

# Data import and preparation

``` r
hundred_sw_clean <- read_excel(
  here("data", "raw", "All Treatments", "combined_raw.xlsx")
) |>
  clean_names() |>
  rename(weed_trt = treatment) |>
  mutate(
    year      = factor(year),
    location  = factor(location),
    site_year = factor(interaction(year, location, drop = TRUE)),
    block     = factor(block),
    weed_trt  = recode(
      weed_trt,
      "RNO" = "Rolled, no control",
      "RIM" = "Rolled, mowing",
      "RIC" = "Rolled, high-residue cultivation",
      "TIM" = "Tilled, mowing",
      "TIC" = "Tilled, cultivation"
    ),
    weed_trt = factor(weed_trt, levels = mow_levels)
  ) |>
  # keep only rows with non-missing 100-seed weight
  filter(!is.na(seed_weight))
  # no unit conversions here; use hundred_seed_weight as recorded


# Quick check
kable(
  head(hundred_sw_clean),
  caption = "All site-years, cleaned (100-seed weight)"
)
```

| id | location | year | weed_trt | block | plot | bean_emergence | bean_biomass | inrow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | seed_weight | site_year |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| CU_B1_P101 | field v | 2023 | Tilled, mowing | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 17.1200 | 2023.field v |
| CU_B1_P102 | field v | 2023 | Tilled, cultivation | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 17.4750 | 2023.field v |
| CU_B1_P103 | field v | 2023 | Rolled, mowing | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 7.840 | 37.5 | 449.93 | 16.7525 | 2023.field v |
| CU_B1_P104 | field v | 2023 | Rolled, no control | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35.0 | 412.59 | 16.1450 | 2023.field v |
| CU_B1_P105 | field v | 2023 | Rolled, high-residue cultivation | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39.0 | 473.79 | 17.0475 | 2023.field v |
| CU_B1_P201 | field v | 2023 | Rolled, high-residue cultivation | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 17.1500 | 2023.field v |

All site-years, cleaned (100-seed weight)

# Model testing and diagnostics

## Exploratory checks

``` r
## Exploratory checks: 100-seed weight (g) by site-year × treatment -------

## 0) Directories for 100-seed-weight tables & figures -------------------
tab_dir_100sw <- here("analysis", "tables", "100-seed-weight")
dir.create(tab_dir_100sw, showWarnings = FALSE, recursive = TRUE)

fig_dir_100sw <- here("analysis", "figs", "100-seed-weight")
dir.create(fig_dir_100sw, showWarnings = FALSE, recursive = TRUE)


## 1) Summary table: 100-seed weight by site-year × treatment ------------

sw_sy_trt_summary <- hundred_sw_clean |>
  mutate(
    site_year = as.factor(site_year),
    weed_trt  = factor(weed_trt, levels = mow_levels)
  ) |>
  group_by(site_year, weed_trt) |>
  summarise(
    n      = sum(!is.na(seed_weight)),
    mean   = mean(seed_weight, na.rm = TRUE),
    median = median(seed_weight, na.rm = TRUE),
    sd     = sd(seed_weight, na.rm = TRUE),
    se     = sd / sqrt(n),
    .groups = "drop"
  ) |>
  arrange(site_year, weed_trt)

# Save summary table
readr::write_csv(
  sw_sy_trt_summary,
  file.path(tab_dir_100sw, "tab_100-seed-weight_site-year_treatment_summary_raw_g.csv")
)

sw_sy_trt_summary |>
  kable(
    digits  = 2,
    caption = "100-seed weight (g) by site-year × treatment (raw means ± SD and SE)"
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g) by site-year × treatment (raw means ± SD and SE)
</caption>

<thead>

<tr>

<th style="text-align:left;">

site_year
</th>

<th style="text-align:left;">

weed_trt
</th>

<th style="text-align:right;">

n
</th>

<th style="text-align:right;">

mean
</th>

<th style="text-align:right;">

median
</th>

<th style="text-align:right;">

sd
</th>

<th style="text-align:right;">

se
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

17.69
</td>

<td style="text-align:right;">

17.69
</td>

<td style="text-align:right;">

0.05
</td>

<td style="text-align:right;">

0.02
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

17.88
</td>

<td style="text-align:right;">

17.88
</td>

<td style="text-align:right;">

0.24
</td>

<td style="text-align:right;">

0.12
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.06
</td>

<td style="text-align:right;">

18.06
</td>

<td style="text-align:right;">

0.21
</td>

<td style="text-align:right;">

0.11
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.63
</td>

<td style="text-align:right;">

18.46
</td>

<td style="text-align:right;">

0.61
</td>

<td style="text-align:right;">

0.31
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.66
</td>

<td style="text-align:right;">

18.95
</td>

<td style="text-align:right;">

0.68
</td>

<td style="text-align:right;">

0.34
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.24
</td>

<td style="text-align:right;">

18.16
</td>

<td style="text-align:right;">

0.33
</td>

<td style="text-align:right;">

0.16
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.27
</td>

<td style="text-align:right;">

18.07
</td>

<td style="text-align:right;">

0.50
</td>

<td style="text-align:right;">

0.25
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.30
</td>

<td style="text-align:right;">

18.20
</td>

<td style="text-align:right;">

0.53
</td>

<td style="text-align:right;">

0.26
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

18.28
</td>

<td style="text-align:right;">

18.45
</td>

<td style="text-align:right;">

0.85
</td>

<td style="text-align:right;">

0.43
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

19.14
</td>

<td style="text-align:right;">

19.12
</td>

<td style="text-align:right;">

0.23
</td>

<td style="text-align:right;">

0.12
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

16.38
</td>

<td style="text-align:right;">

16.34
</td>

<td style="text-align:right;">

0.23
</td>

<td style="text-align:right;">

0.12
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

16.80
</td>

<td style="text-align:right;">

16.74
</td>

<td style="text-align:right;">

0.29
</td>

<td style="text-align:right;">

0.15
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

17.50
</td>

<td style="text-align:right;">

17.42
</td>

<td style="text-align:right;">

0.49
</td>

<td style="text-align:right;">

0.25
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

17.21
</td>

<td style="text-align:right;">

17.09
</td>

<td style="text-align:right;">

0.45
</td>

<td style="text-align:right;">

0.23
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

17.25
</td>

<td style="text-align:right;">

17.32
</td>

<td style="text-align:right;">

0.25
</td>

<td style="text-align:right;">

0.13
</td>

</tr>

</tbody>

</table>

``` r
## 2) Faceted boxplot: all site-years -----------------------------------

fig_100sw_box_sy <- hundred_sw_clean |>
  mutate(
    weed_trt  = factor(weed_trt, levels = mow_levels),
    site_year = as.factor(site_year)
  ) |>
  ggplot(aes(x = weed_trt, y = seed_weight, fill = weed_trt)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "black") +
  geom_jitter(width = 0.12, height = 0, alpha = 0.4, size = 1.8, color = "grey30") +
  facet_wrap(~ site_year, nrow = 1) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_comma_cult) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x     = NULL,
    y     = "100-seed weight (g)",
    title = "100-seed weight by treatment across site-years"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text  = element_text(face = "bold")
  )

fig_100sw_box_sy  # print in the Rmd
```

![](analysis/fig/bean_100_weight-unnamed-chunk-3-1.png)<!-- -->

``` r
# Save figure
ggsave(
  filename = file.path(fig_dir_100sw, "fig_100-seed-weight_box_by-site-year_raw_g.png"),
  plot     = fig_100sw_box_sy,
  width    = 12,
  height   = 5.5,
  dpi      = 300
)
```

## Model selection

``` r
### Model testing / selection for 100-seed weight (g) --------------------

options(contrasts = c("contr.sum", "contr.poly"))

# 1) Fit Gaussian LMMs with REML (for estimation) ------------------------

sw_int <- lmer(
  seed_weight ~ weed_trt * site_year + (1 | site_year:block),
  data = hundred_sw_clean
)

sw_add <- lmer(
  seed_weight ~ weed_trt + site_year + (1 | site_year:block),
  data = hundred_sw_clean
)

# 2) Likelihood-ratio test and ML-based AIC via anova() ------------------
# anova() refits models with REML = FALSE internally for the comparison

lrt_sw <- anova(sw_add, sw_int)  # ML fits used here

# Extract ML AICs (these are the ones you see in the printed table)
AIC_add_sw   <- lrt_sw$AIC[1]
AIC_int_sw   <- lrt_sw$AIC[2]
deltaAIC_sw  <- AIC_add_sw - AIC_int_sw   # > 0 => interaction has lower AIC

# Extract LRT p-value
p_int_sw <- lrt_sw$`Pr(>Chisq)`[2]

# 3) Classify evidence for interaction ----------------------------------

p_strong_sw    <- 0.01
p_none_sw      <- 0.20
dAIC_strong_sw <- 4

interaction_class_sw <- dplyr::case_when(
  p_int_sw < p_strong_sw | deltaAIC_sw >= dAIC_strong_sw ~ "interaction",
  p_int_sw > p_none_sw  & abs(deltaAIC_sw) < 2           ~ "additive",
  TRUE                                                   ~ "gray_zone"
)

primary_model_name_sw <- dplyr::case_when(
  interaction_class_sw == "interaction" ~ "Interaction: weed_trt * site_year",
  TRUE                                  ~ "Additive: weed_trt + site_year"
)

# 4) Final model used for emmeans / plots (REML fits) --------------------

sw.lmer <- if (primary_model_name_sw == "Interaction: weed_trt * site_year") {
  sw_int
} else {
  sw_add
}

family_structure_sw <- "Gaussian LMM (identity link)"

# 5) AIC table for reporting (using ML AICs from lrt_sw) -----------------

aic_sw_out <- tibble(
  model = c(
    "Additive: weed_trt + site_year",
    "Interaction: weed_trt * site_year"
  ),
  AIC = c(AIC_add_sw, AIC_int_sw)
) |>
  dplyr::mutate(
    deltaAIC                 = AIC - min(AIC),
    Selected                 = dplyr::if_else(model == primary_model_name_sw, "Yes", ""),
    Evidence_for_interaction = interaction_class_sw
  )

kable(
  aic_sw_out,
  digits  = 2,
  caption = "100-seed weight (g): additive vs interaction (Gaussian LMM, ML AIC)"
) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g): additive vs interaction (Gaussian LMM, ML AIC)
</caption>

<thead>

<tr>

<th style="text-align:left;">

model
</th>

<th style="text-align:right;">

AIC
</th>

<th style="text-align:right;">

deltaAIC
</th>

<th style="text-align:left;">

Selected
</th>

<th style="text-align:left;">

Evidence_for_interaction
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Additive: weed_trt + site_year
</td>

<td style="text-align:right;">

90.05
</td>

<td style="text-align:right;">

0.00
</td>

<td style="text-align:left;">

Yes
</td>

<td style="text-align:left;">

gray_zone
</td>

</tr>

<tr>

<td style="text-align:left;">

Interaction: weed_trt \* site_year
</td>

<td style="text-align:right;">

90.72
</td>

<td style="text-align:right;">

0.67
</td>

<td style="text-align:left;">

</td>

<td style="text-align:left;">

gray_zone
</td>

</tr>

</tbody>

</table>

## Diagnostics (residuals, dispersion, zero inflation)

``` r
# 1) Quick reminder + diagnostics ---------------------------------------

cat(
  "\nSelected primary model for 100-seed weight (used in emmeans/plots):\n  ",
  primary_model_name_sw,
  sprintf(
    "  [LRT p = %.3f; ΔAIC (add - int) = %.2f; class = %s]\n",
    p_int_sw, deltaAIC_sw, interaction_class_sw
  )
)
```

    ## 
    ## Selected primary model for 100-seed weight (used in emmeans/plots):
    ##    Additive: weed_trt + site_year   [LRT p = 0.053; ΔAIC (add - int) = -0.67; class = gray_zone]

``` r
set.seed(123)
res_sw <- DHARMa::simulateResiduals(sw.lmer)
plot(res_sw)
```

![](analysis/fig/bean_100_weight-unnamed-chunk-5-1.png)<!-- -->

``` r
DHARMa::testDispersion(sw.lmer)
```

![](analysis/fig/bean_100_weight-unnamed-chunk-5-2.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.89294, p-value = 0.552
    ## alternative hypothesis: two.sided

``` r
car::Anova(sw.lmer, type = 3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: seed_weight
    ##                 Chisq Df Pr(>Chisq)    
    ## (Intercept) 87159.654  1  < 2.2e-16 ***
    ## weed_trt       27.227  4  1.789e-05 ***
    ## site_year     103.559  2  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Summary tables

## Post-hoc treatment means and CLDs

``` r
### 100-seed weight (g) with Fisher's LSD CLDs

# Estimated marginal means for weed_trt
emm_sw <- emmeans(sw.lmer, ~ weed_trt)

# Tidy emmeans (adds ci_low / ci_high and enforces treatment order)
emm_sw_df <- tidy_emm(emm_sw, ref_levels = mow_levels) |>
  as_tibble()

# Compact letter display (Fisher's LSD, no adjustment; "a" = highest)
cld_sw <- cld(
  emm_sw,
  adjust  = "none",
  Letters = letters,
  sort    = TRUE,   # explicit
  reversed = TRUE   # "a" = highest mean
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# Join emmeans + CLDs and format for reporting
emm_sw_df |>
  left_join(cld_sw, by = "weed_trt") |>
  select(weed_trt, emmean, SE, ci_low, ci_high, .group) |>
  mutate(across(c(emmean, SE, ci_low, ci_high), ~ round(.x, 2))) |>
  kable(
    caption   = "Estimated 100-seed weight (g) with 95% CI and Fisher's LSD group letters",
    col.names = c("Treatment", "Mean", "SE", "Lower CI", "Upper CI", "Group")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Estimated 100-seed weight (g) with 95% CI and Fisher’s LSD group letters
</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment
</th>

<th style="text-align:right;">

Mean
</th>

<th style="text-align:right;">

SE
</th>

<th style="text-align:right;">

Lower CI
</th>

<th style="text-align:right;">

Upper CI
</th>

<th style="text-align:left;">

Group
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

17.44
</td>

<td style="text-align:right;">

0.14
</td>

<td style="text-align:right;">

17.17
</td>

<td style="text-align:right;">

17.71
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

17.65
</td>

<td style="text-align:right;">

0.14
</td>

<td style="text-align:right;">

17.38
</td>

<td style="text-align:right;">

17.92
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:right;">

0.14
</td>

<td style="text-align:right;">

17.68
</td>

<td style="text-align:right;">

18.22
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.04
</td>

<td style="text-align:right;">

0.14
</td>

<td style="text-align:right;">

17.77
</td>

<td style="text-align:right;">

18.31
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.35
</td>

<td style="text-align:right;">

0.14
</td>

<td style="text-align:right;">

18.08
</td>

<td style="text-align:right;">

18.63
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

## Global response summary (ANOVA, AIC, model info)

``` r
#### Global response summary table: 100-seed weight (g)

## 0) Directory for all 100-seed-weight tables ---------------------------

tab_dir_sw <- here("analysis", "tables", "100-seed-weight")
dir.create(tab_dir_sw, showWarnings = FALSE, recursive = TRUE)


# 1) P-value summary (Location, Treatment, Interaction) ------------------

anova_sw <- Anova(sw.lmer, type = 3)

anova_sw_df <- anova_sw |>
  as.data.frame() |>
  tibble::rownames_to_column("Effect")

p_site_sw <- anova_sw_df$`Pr(>Chisq)`[anova_sw_df$Effect == "site_year"]
p_trt_sw  <- anova_sw_df$`Pr(>Chisq)`[anova_sw_df$Effect == "weed_trt"]
p_int_sw3 <- anova_sw_df$`Pr(>Chisq)`[anova_sw_df$Effect == "weed_trt:site_year"]

pvals_sw <- tibble(
  Effect = c(
    "Location (site_year)",
    "Treatment (weed_trt)"
  ),
  p_raw = c(p_site_sw, p_trt_sw)
)

# Only include the interaction row if the *selected* primary model uses the interaction
if (primary_model_name_sw == "Interaction: weed_trt * site_year") {
  pvals_sw <- bind_rows(
    pvals_sw,
    tibble(
      Effect = "Location × Treatment",
      p_raw  = p_int_sw3
    )
  )
}

pvals_sw <- pvals_sw |>
  mutate(
    `P-value` = case_when(
      p_raw < 0.001 ~ "<0.001",
      p_raw < 0.01  ~ "<0.01",
      TRUE          ~ sprintf("%.3f", p_raw)
    )
  ) |>
  select(Effect, `P-value`)

# Save ANOVA p-value summary
readr::write_csv(
  pvals_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_Anova_pvals.csv")
)


## 1b) Likelihood-ratio test: additive vs interaction --------------------

lrt_table_sw <- tibble(
  Test  = "LRT (additive vs interaction)",
  p_raw = p_int_sw
) |>
  mutate(
    `P-value` = case_when(
      p_raw < 0.001 ~ "<0.001",
      p_raw < 0.01  ~ "<0.01",
      TRUE          ~ sprintf("%.3f", p_raw)
    )
  ) |>
  select(Test, `P-value`)

# Save LRT summary
readr::write_csv(
  lrt_table_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_LRT_add-vs-int.csv")
)

lrt_table_sw |>
  kable(
    caption   = "100-seed weight (g): Likelihood-ratio test comparing additive vs interaction models",
    col.names = c("Test", "P-value")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g): Likelihood-ratio test comparing additive vs
interaction models
</caption>

<thead>

<tr>

<th style="text-align:left;">

Test
</th>

<th style="text-align:left;">

P-value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

LRT (additive vs interaction)
</td>

<td style="text-align:left;">

0.053
</td>

</tr>

</tbody>

</table>

``` r
## 2) Location block: site-year means (model + raw) ----------------------

emm_loc_sw <- emmeans(
  sw.lmer,
  ~ site_year
)

emm_loc_sw_df <- as_tibble(emm_loc_sw) |>
  mutate(
    site_year  = as.factor(site_year),
    model_mean = emmean   # g per 100 seeds
  ) |>
  select(site_year, model_mean)

cld_loc_sw <- cld(
  emm_loc_sw,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    site_year = as.factor(site_year),
    loc_CLD   = str_trim(.group)
  ) |>
  select(site_year, loc_CLD)

raw_loc_sw <- hundred_sw_clean |>
  group_by(site_year) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(site_year = as.factor(site_year))

loc_summary_sw <- emm_loc_sw_df |>
  left_join(cld_loc_sw, by = "site_year") |>
  left_join(raw_loc_sw, by = "site_year") |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = loc_CLD
  ) |>
  arrange(site_year)

# Save location summary
readr::write_csv(
  loc_summary_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_location_means_CLD.csv")
)

loc_summary_sw |>
  kable(
    caption   = "100-seed weight (g): location (site-year) means with CLDs",
    col.names = c("Site-year", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g): location (site-year) means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Site-year
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:right;">

18.19
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.19
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:right;">

18.45
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.45
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:right;">

17.03
</td>

<td style="text-align:left;">

b
</td>

<td style="text-align:right;">

17.03
</td>

<td style="text-align:left;">

b
</td>

</tr>

</tbody>

</table>

``` r
## 3) Treatment block: means (model + raw) -------------------------------

emm_sw_trt <- emmeans(
  sw.lmer,
  ~ weed_trt
)

emm_trt_sw_df <- as_tibble(emm_sw_trt) |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    model_mean = emmean
  ) |>
  select(weed_trt, model_mean)

cld_sw_trt <- cld(
  emm_sw_trt,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    trt_CLD  = str_trim(.group)
  ) |>
  select(weed_trt, trt_CLD)

raw_trt_sw <- hundred_sw_clean |>
  group_by(weed_trt) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(weed_trt = factor(weed_trt, levels = mow_levels))

trt_summary_sw <- emm_trt_sw_df |>
  left_join(cld_sw_trt,  by = "weed_trt") |>
  left_join(raw_trt_sw,  by = "weed_trt") |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = trt_CLD
  ) |>
  arrange(weed_trt)

# Save treatment summary
readr::write_csv(
  trt_summary_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_treatment_means_CLD.csv")
)

trt_summary_sw |>
  kable(
    caption   = "100-seed weight (g): treatment means with CLDs",
    col.names = c("Treatment", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g): treatment means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

17.44
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

17.44
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

17.65
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

17.65
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.04
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.04
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.35
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.35
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

``` r
## 4) Interaction block: site-year × treatment means ---------------------

emm_sy_sw <- emmeans(
  sw.lmer,
  ~ weed_trt | site_year
)

emm_sy_sw_df <- as_tibble(emm_sy_sw) |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    site_year  = as.factor(site_year),
    model_mean = emmean
  ) |>
  select(site_year, weed_trt, model_mean)

cld_sy_sw <- cld(
  emm_sy_sw,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    weed_trt  = factor(weed_trt, levels = mow_levels),
    site_year = as.factor(site_year),
    int_CLD   = str_trim(.group)
  ) |>
  select(site_year, weed_trt, int_CLD)

raw_sy_sw <- hundred_sw_clean |>
  group_by(site_year, weed_trt) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    site_year = as.factor(site_year),
    weed_trt  = factor(weed_trt, levels = mow_levels)
  )

int_summary_sw <- emm_sy_sw_df |>
  left_join(cld_sy_sw, by = c("site_year", "weed_trt")) |>
  left_join(raw_sy_sw, by = c("site_year", "weed_trt")) |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = int_CLD
  ) |>
  arrange(site_year, weed_trt)

# Save interaction summary
readr::write_csv(
  int_summary_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_site-year_treatment_means_CLD.csv")
)

int_summary_sw |>
  kable(
    caption   = "100-seed weight (g): site-year × treatment means with CLDs",
    col.names = c(
      "Site-year", "Treatment",
      "Model mean", "Model CLD",
      "Raw mean",   "Raw CLD"
    )
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

100-seed weight (g): site-year × treatment means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Site-year
</th>

<th style="text-align:left;">

Treatment
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

17.74
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

17.69
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

17.88
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

18.25
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

18.06
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.34
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.63
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.65
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.66
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

18.00
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

18.24
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

18.21
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

18.27
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

18.51
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

18.30
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.60
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.28
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.91
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

19.14
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

16.58
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

16.38
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

16.79
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

16.80
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

17.09
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

17.50
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

17.18
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

17.21
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

17.50
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

17.25
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

``` r
## 5) Model-info row for global response summary -------------------------

model_info_sw <- tibble::tibble(
  response_label    = "100-seed weight (g), model-predicted means*",
  family_structure  = family_structure_sw,
  fixed_effects     = primary_model_name_sw,
  random_effects    = "(1 | site_year:block)",
  AIC_additive      = round(AIC_add_sw, 2),
  AIC_interaction   = round(AIC_int_sw, 2),
  deltaAIC_add_int  = round(deltaAIC_sw, 2),
  LRT_p_int_raw     = p_int_sw,
  LRT_p_int_label   = dplyr::case_when(
    p_int_sw < 0.001 ~ "<0.001",
    p_int_sw < 0.01  ~ "<0.01",
    TRUE             ~ sprintf("%.3f", p_int_sw)
  ),
  interaction_class = interaction_class_sw
)

readr::write_csv(
  model_info_sw,
  file.path(tab_dir_sw, "tab_100-seed-weight_model-info.csv")
)


## 6) OPTIONAL: Combine all 100-seed-weight tables into one Excel workbook
## (requires: library(writexl) in your Packages chunk)

sw_tables <- list(
  Anova_pvals            = pvals_sw,
  LRT_add_vs_int         = lrt_table_sw,
  Location_means_CLD     = loc_summary_sw,
  Treatment_means_CLD    = trt_summary_sw,
  SiteYear_trt_means     = int_summary_sw,
  Model_info             = model_info_sw
)

write_xlsx(
  sw_tables,
  path = file.path(tab_dir_sw, "100-seed-weight_all-tables.xlsx")
)
```

``` r
## 1) P-value summary (Location, Treatment, Interaction) -----------------

# Type-III tests for additive model (weed_trt + site_year)
anova_sw <- Anova(sw.lmer, type = 3)

anova_sw_df <- anova_sw |>
  as.data.frame() |>
  tibble::rownames_to_column("Effect")

# LRT for interaction (additive vs interaction models)
anova_interaction_sw <- anova(sw_add, sw_int)

pvals_sw <- tibble(
  Effect = c(
    "Location (site_year)",
    "Treatment (weed_trt)",
    "Location × Treatment (LRT)"
  ),
  p_raw  = c(
    anova_sw_df$`Pr(>Chisq)`[anova_sw_df$Effect == "site_year"],
    anova_sw_df$`Pr(>Chisq)`[anova_sw_df$Effect == "weed_trt"],
    anova_interaction_sw$`Pr(>Chisq)`[2]   # LRT p-value
  )
) |>
  mutate(
    `P-value` = case_when(
      p_raw < 0.001 ~ "<0.001",
      p_raw < 0.01  ~ "<0.01",
      TRUE          ~ sprintf("%.3f", p_raw)
    )
  ) |>
  select(Effect, `P-value`)

pvals_sw |>
  kable(
    caption   = "Soybean 100-seed weight (g): P-values for location, treatment, and location × treatment (interaction p from LRT).",
    col.names = c("Effect", "P-value")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Soybean 100-seed weight (g): P-values for location, treatment, and
location × treatment (interaction p from LRT).
</caption>

<thead>

<tr>

<th style="text-align:left;">

Effect
</th>

<th style="text-align:left;">

P-value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Location (site_year)
</td>

<td style="text-align:left;">

\<0.001
</td>

</tr>

<tr>

<td style="text-align:left;">

Treatment (weed_trt)
</td>

<td style="text-align:left;">

\<0.001
</td>

</tr>

<tr>

<td style="text-align:left;">

Location × Treatment (LRT)
</td>

<td style="text-align:left;">

0.053
</td>

</tr>

</tbody>

</table>

``` r
## 2) Location block: site-year means (model + raw) ----------------------

# Model-based emmeans by site_year
emm_loc_sw <- emmeans(sw.lmer, ~ site_year)

emm_loc_sw_df <- tidy_emm(emm_loc_sw) |>
  as_tibble() |>
  mutate(
    site_year  = as.factor(site_year),
    model_mean = emmean
  ) |>
  select(site_year, model_mean)

# CLDs for site_year (a = highest)
cld_loc_sw <- cld(
  emm_loc_sw,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    site_year = as.factor(site_year),
    loc_CLD   = str_trim(.group)
  ) |>
  select(site_year, loc_CLD)

# Raw means by site_year
raw_loc_sw <- hundred_sw_clean |>
  group_by(site_year) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(site_year = as.factor(site_year))

loc_summary_sw <- emm_loc_sw_df |>
  left_join(cld_loc_sw, by = "site_year") |>
  left_join(raw_loc_sw, by = "site_year") |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = loc_CLD  # use same letters for model + raw
  ) |>
  arrange(site_year)

loc_summary_sw |>
  kable(
    caption   = "Soybean 100-seed weight (g): location (site-year) means with CLDs",
    col.names = c("Site-year", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Soybean 100-seed weight (g): location (site-year) means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Site-year
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:right;">

18.19
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.19
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:right;">

18.45
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.45
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:right;">

17.03
</td>

<td style="text-align:left;">

b
</td>

<td style="text-align:right;">

17.03
</td>

<td style="text-align:left;">

b
</td>

</tr>

</tbody>

</table>

``` r
## 3) Treatment block: means (model + raw) -------------------------------

# emmeans for treatment
emm_sw <- emmeans(sw.lmer, ~ weed_trt)

emm_trt_sw_df <- tidy_emm(emm_sw, ref_levels = mow_levels) |>
  as_tibble() |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    model_mean = emmean
  ) |>
  select(weed_trt, model_mean)

# CLDs for treatment (a = highest)
cld_sw <- cld(
  emm_sw,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    trt_CLD  = str_trim(.group)
  ) |>
  select(weed_trt, trt_CLD)

# Raw means by treatment
raw_trt_sw <- hundred_sw_clean |>
  group_by(weed_trt) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(weed_trt = factor(weed_trt, levels = mow_levels))

trt_summary_sw <- emm_trt_sw_df |>
  left_join(cld_sw, by = "weed_trt") |>
  left_join(raw_trt_sw, by = "weed_trt") |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = trt_CLD
  ) |>
  arrange(weed_trt)

trt_summary_sw |>
  kable(
    caption   = "Soybean 100-seed weight (g): treatment means with CLDs",
    col.names = c("Treatment", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Soybean 100-seed weight (g): treatment means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Treatment
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

17.44
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

17.44
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

17.65
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

17.65
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.04
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.04
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.35
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.35
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

``` r
## 4) Interaction block: site-year × treatment means ---------------------

# Model emmeans by treatment within site_year
emm_sy_sw <- emmeans(sw.lmer, ~ weed_trt | site_year)

emm_sy_sw_df <- tidy_emm(emm_sy_sw, ref_levels = mow_levels) |>
  as_tibble() |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    site_year  = as.factor(site_year),
    model_mean = emmean
  ) |>
  select(site_year, weed_trt, model_mean)

# CLDs within each site_year (a = highest within that site_year)
cld_sy_sw <- cld(
  emm_sy_sw,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    weed_trt  = factor(weed_trt, levels = mow_levels),
    site_year = as.factor(site_year),
    int_CLD   = str_trim(.group)
  ) |>
  select(site_year, weed_trt, int_CLD)

# Raw means by site_year × treatment
raw_sy_sw <- hundred_sw_clean |>
  group_by(site_year, weed_trt) |>
  summarise(
    raw_mean = mean(seed_weight, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    site_year = as.factor(site_year),
    weed_trt  = factor(weed_trt, levels = mow_levels)
  )

int_summary_sw <- emm_sy_sw_df |>
  left_join(cld_sy_sw, by = c("site_year", "weed_trt")) |>
  left_join(raw_sy_sw, by = c("site_year", "weed_trt")) |>
  mutate(
    model_mean = round(model_mean, 2),
    raw_mean   = round(raw_mean, 2),
    raw_CLD    = int_CLD
  ) |>
  arrange(site_year, weed_trt)

int_summary_sw |>
  kable(
    caption   = "Soybean 100-seed weight (g): site-year × treatment means with CLDs",
    col.names = c(
      "Site-year", "Treatment",
      "Model mean", "Model CLD",
      "Raw mean",   "Raw CLD"
    )
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Soybean 100-seed weight (g): site-year × treatment means with CLDs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Site-year
</th>

<th style="text-align:left;">

Treatment
</th>

<th style="text-align:right;">

Model mean
</th>

<th style="text-align:left;">

Model CLD
</th>

<th style="text-align:right;">

Raw mean
</th>

<th style="text-align:left;">

Raw CLD
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

17.74
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

17.69
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

17.95
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

17.88
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

18.25
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

18.06
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.34
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.63
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.65
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

18.66
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

18.00
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

18.24
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

18.21
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

18.27
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

18.51
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

18.30
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

18.60
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

18.28
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

18.91
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

19.14
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, no control
</td>

<td style="text-align:right;">

16.58
</td>

<td style="text-align:left;">

d
</td>

<td style="text-align:right;">

16.38
</td>

<td style="text-align:left;">

d
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, mowing
</td>

<td style="text-align:right;">

16.79
</td>

<td style="text-align:left;">

cd
</td>

<td style="text-align:right;">

16.80
</td>

<td style="text-align:left;">

cd
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled, high-residue cultivation
</td>

<td style="text-align:right;">

17.09
</td>

<td style="text-align:left;">

bc
</td>

<td style="text-align:right;">

17.50
</td>

<td style="text-align:left;">

bc
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, mowing
</td>

<td style="text-align:right;">

17.18
</td>

<td style="text-align:left;">

ab
</td>

<td style="text-align:right;">

17.21
</td>

<td style="text-align:left;">

ab
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled, cultivation
</td>

<td style="text-align:right;">

17.50
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

17.25
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

# Figures

## Model-predicted means (g 100 seeds )

``` r
## 0) Directory for 100-seed-weight figures ------------------------------
fig_dir_sw <- here("analysis", "figs", "100-seed-weight")
dir.create(fig_dir_sw, showWarnings = FALSE, recursive = TRUE)

## Figure: Soybean 100-seed weight by weed management treatment
## (model-predicted means ± SE, model CLDs)

# 1) Model-based emmeans by treatment -----------------------------------
emm_sw_trt <- emmeans(
  sw.lmer,
  ~ weed_trt
)

# Tidy emmeans and prepare SE-based error bars ---------------------------
emm_sw_trt_df <- as_tibble(emm_sw_trt) |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    mean     = emmean,            # LMM-predicted marginal mean (g per 100 seeds)
    ymin     = pmax(mean - SE, 0),
    ymax     = mean + SE
  )

# 2) Model-based CLDs for treatment main effect --------------------------
cld_sw_trt <- cld(
  emm_sw_trt,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE   # "a" = highest 100-seed weight group(s)
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# 3) Join model-predicted means + CLDs ----------------------------------
plot_df_sw_model <- emm_sw_trt_df |>
  left_join(cld_sw_trt, by = "weed_trt")

# 4) Plot ---------------------------------------------------------------
fig_sw_model <- ggplot(
  plot_df_sw_model,
  aes(x = weed_trt, y = mean, fill = weed_trt)
) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust    = 0,
    fontface = "bold",
    size     = 6
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_comma_cult) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x       = NULL,
    y       = "Soybean 100-seed weight (g)",
    title   = "Soybean 100-seed weight by weed management",
    caption = "Bars show LMM-predicted marginal means (g) ± SE; letters denote Fisher’s LSD groups (α = 0.05)."
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title   = element_text(face = "bold"),
    plot.caption = element_text(size = 9, hjust = 0)
  )

fig_sw_model  # print in the Rmd
```

![](analysis/fig/bean_100_weight-unnamed-chunk-9-1.png)<!-- -->

``` r
# 5) Save figure ---------------------------------------------------------
ggsave(
  filename = file.path(fig_dir_sw, "fig_100-seed-weight_total_model_g.png"),
  plot     = fig_sw_model,
  width    = 9,
  height   = 5.5,
  dpi      = 300
)
```

## Raw means (g 100 seeds )

``` r
## Raw means (g) ± SE (letters = model-based CLDs) -----------------------

## 0) Directory for 100-seed-weight figures ------------------------------
fig_dir_sw <- here("analysis", "figs", "100-seed-weight")
dir.create(fig_dir_sw, showWarnings = FALSE, recursive = TRUE)

## Figure: Soybean 100-seed weight by weed management treatment
## (raw means ± SE, model-based CLDs)

# 1) Raw means and SE by treatment --------------------------------------
raw_sw_summary <- hundred_sw_clean |>
  group_by(weed_trt) |>
  summarise(
    n    = n(),
    mean = mean(seed_weight, na.rm = TRUE),
    sd   = sd(seed_weight, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    ymin     = pmax(mean - se, 0),
    ymax     = mean + se
  )

# 2) Model-based CLDs for treatment main effect --------------------------
emm_sw_trt <- emmeans(
  sw.lmer,
  ~ weed_trt
)

cld_sw_trt <- cld(
  emm_sw_trt,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE   # "a" = highest 100-seed weight group(s)
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# 3) Join raw means + model CLDs ----------------------------------------
plot_df_sw_raw <- raw_sw_summary |>
  left_join(cld_sw_trt, by = "weed_trt")

# 4) Plot ----------------------------------------------------------------
fig_sw_raw <- ggplot(
  plot_df_sw_raw,
  aes(x = weed_trt, y = mean, fill = weed_trt)
) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust    = 0,
    fontface = "bold",
    size     = 6
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_comma_cult) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x       = NULL,
    y       = "Soybean 100-seed weight (g)",
    title   = "Soybean 100-seed weight by weed management",
    caption = "Bars show raw means (g) ± SE; letters denote model-based Fisher’s LSD groups (α = 0.05)."
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title   = element_text(face = "bold"),
    plot.caption = element_text(size = 9, hjust = 0)
  )

fig_sw_raw  # print in the Rmd
```

![](analysis/fig/bean_100_weight-unnamed-chunk-10-1.png)<!-- -->

``` r
# 5) Save figure ---------------------------------------------------------
ggsave(
  filename = file.path(fig_dir_sw, "fig_100-seed-weight_total_raw_g.png"),
  plot     = fig_sw_raw,
  width    = 9,
  height   = 5.5,
  dpi      = 300
)
```
