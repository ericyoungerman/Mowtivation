Soybean crop biomass
================

- [Question and experimental design](#question-and-experimental-design)
- [Analysis overview](#analysis-overview)
- [Setup](#setup)
- [Import and validate](#import-and-validate)
- [Prepare subplot observations](#prepare-subplot-observations)
- [Describe the observations](#describe-the-observations)
- [Fit the primary model](#fit-the-primary-model)
- [Evaluate random-effects
  structure](#evaluate-random-effects-structure)
- [Sensitivity to the random-effects
  structure](#sensitivity-to-the-random-effects-structure)
- [Check model assumptions](#check-model-assumptions)
- [Estimated marginal means](#estimated-marginal-means)
- [Compare management systems](#compare-management-systems)
- [Fisher’s LSD and letter groups](#fishers-lsd-and-letter-groups)
- [Model-estimated biomass figures](#model-estimated-biomass-figures)
- [Alternative figures: observed means and descriptive
  SEs](#alternative-figures-observed-means-and-descriptive-ses)
- [Site-year comparisons](#site-year-comparisons)
- [Results summary](#results-summary)
- [Appendix: sensitivity to residual
  variance](#appendix-sensitivity-to-residual-variance)
- [Reproducibility and remaining
  limitations](#reproducibility-and-remaining-limitations)

# Question and experimental design

Did soybean aboveground biomass differ among five management systems,
and did system differences vary among experiments (site-years)?

The experiment evaluated five management systems in a split-plot design
across three site-years, with four blocks per site-year (60 subplot
observations). Cereal rye termination method was assigned to whole plots
within blocks: moldboard plowing (tilled) or roller-crimping (rolled).
Within tilled whole plots, subplots received cultivation (TIC) or mowing
(TIM). Within rolled whole plots, subplots received no supplemental weed
control (RNO), mowing (RIM), or high-residue cultivation (RIC). One
soybean biomass observation was collected per subplot. Thus, subplots
are the observational units, while comparisons between systems with
different termination methods depend on whole-plot replication. The five
systems are not a fully crossed factorial treatment.

In 2023, an operations error caused one subplot randomized to RNO to
receive RIC (high-residue cultivation). All other treatments were
implemented as assigned. This error did not alter the termination whole
plots: TIM and TIC remained together, and the three rolled-system
subplots remained within their roller-crimped whole plot. The primary
analysis retains the original randomized treatment labels rather than
reclassifying the affected subplot. The processed data retain one record
for each randomized treatment per block; the affected subplot’s identity
is not independently verified here against a randomization map. RIC
received only one late cultivation event rather than the intended
repeated cultivations, so interpretation concerns the management as
implemented.

# Analysis overview

This analysis evaluates aboveground soybean crop biomass (kg/ha) across
three site-years. Management system and its interaction with site-year
are included as fixed effects. The random-effects structure reflects the
split-plot randomization, with termination methods assigned to whole
plots and management systems evaluated at the subplot level. Comparisons
between termination methods must therefore account for whole-plot
replication.

Observed means and variability describe the field measurements.
Model-based estimated marginal means, standard errors, and comparisons
provide inference while accounting for the experimental design. Figures
using both approaches are retained for transparent reporting.

# Setup

``` r
library(tidyverse)
library(here)
library(conflicted)
library(WrensBookshelf)
library(lmerTest)
library(emmeans)
# multcomp supplies cld(); multcompView supplies its letter-display algorithm.
# nlme is used via explicit namespace calls in the variance sensitivity check.
for (pkg in c("multcomp", "multcompView", "nlme")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

method_levels <- c("Roller-crimper", "Moldboard plow")
treatment_levels <- c("RNO", "RIM", "RIC", "TIM", "TIC")
treatment_labels <- c(
  RNO = "Rolled,\nno control",
  RIM = "Rolled,\nmowing",
  RIC = "Rolled,\nhigh-residue cultivation",
  TIM = "Tilled,\nmowing",
  TIC = "Tilled,\ncultivation"
)

experiment_labels <- c(
  V_2023 = "2023 – Field V",
  O2E_2024 = "2024 – Field O2 East",
  O2W_2024 = "2024 – Field O2 West"
)

experiment_colors <- setNames(
  WrensBookshelf::WB_brewer(
    name = "WhatWellBuild", n = 5, type = "discrete"
  )[c(1, 3, 5)],
  names(experiment_labels)
)
method_lookup <- c(
  RNO = "Roller-crimper", RIM = "Roller-crimper",
  RIC = "Roller-crimper", TIM = "Moldboard plow",
  TIC = "Moldboard plow"
)

treatment_colors <- setNames(
  WrensBookshelf::WB_brewer(name = "WhatWellBuild", n = 5, type = "discrete"),
  treatment_levels
)

theme_set(theme_classic(base_size = 16))
figure_caption_theme <- theme(
  plot.caption = element_text(size = 10, hjust = 0, lineheight = 1.1)
)
lsd_caption <- paste(
  "Letters: mixed-model Fisher's LSD comparisons at α = 0.05.",
  "Shared letters indicate no significant difference detected, not statistical equivalence.",
  sep = "\n"
)

# Plain tables render in both HTML and GitHub Markdown.
show_table <- function(x, caption, digits = 1) {
  knitr::kable(x, digits = digits, caption = caption)
}

# Use the same degrees-of-freedom method in ANOVA and estimated means.
emm_options(lmer.df = "satterthwaite")

fmt_biomass <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
fmt_p <- function(x) {
  if (x < 0.001) "< 0.001" else
    paste("=", formatC(x, format = "f", digits = 3))
}

tab_dir <- here("output", "tables", "soybean_crop_biomass")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
```

# Import and validate

The processed dataset is generated separately in `data_prep.Rmd`; this
analysis uses its biomass values in kg/ha without additional
conversions. Checks require one record per treatment in each of four
blocks within each site-year and ensure that the response is finite and
nonnegative.

``` r
processed_file <- here("data", "processed", "combined_processed.rds")
processed_data <- readRDS(processed_file)
# Identify the exact input without repeating preparation or unit conversions.
show_table(
  tibble(input = "data/processed/combined_processed.rds",
         MD5 = unname(tools::md5sum(processed_file))),
  "Processed-data provenance"
)
```

| input                                 | MD5                              |
|:--------------------------------------|:---------------------------------|
| data/processed/combined_processed.rds | 217fbf0aa2d75c95b9463ee4a0f0bbb8 |

Processed-data provenance

``` r
required <- c(
  "site_year", "block", "treatment", "bean_biomass_kg_ha"
)

missing_columns <- setdiff(required, names(processed_data))

if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

plot_data <- processed_data |>
  select(all_of(required)) |>
  mutate(across(c(site_year, block, treatment), as.character))

if (anyNA(plot_data[c("site_year", "block", "treatment")]) ||
    any(vapply(
      plot_data[c("site_year", "block", "treatment")],
      function(x) any(trimws(x) == ""),
      logical(1)
    ))) {
  stop("Missing or blank site-year, block, or treatment identifiers.")
}

if (!setequal(unique(plot_data$site_year), names(experiment_labels))) {
  stop("Site-year identifiers do not match experiment_labels.")
}

if (!setequal(unique(plot_data$treatment), treatment_levels)) {
  stop("Treatment identifiers do not match treatment_levels.")
}

y <- plot_data$bean_biomass_kg_ha

if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0)) {
  stop("Biomass must contain finite, nonnegative numeric values with no missing data.")
}

plot_counts <- plot_data |>
  count(site_year, block, treatment)

block_counts <- plot_counts |>
  count(site_year, block, name = "n_treatments")

experiment_counts <- block_counts |>
  count(site_year, name = "n_blocks")

if (any(plot_counts$n != 1L) ||
    any(block_counts$n_treatments != 5L) ||
    any(experiment_counts$n_blocks != 4L)) {
  stop(
    "Expected one record per treatment, five treatments per block, ",
    "and four blocks per site-year."
  )
}

plot_data <- plot_data |>
  mutate(
    site_year = factor(site_year, levels = names(experiment_labels)),
    treatment = factor(treatment, levels = treatment_levels),
    primary_method = factor(
      unname(method_lookup[as.character(treatment)]),
      levels = method_levels
    ),
    block_id = interaction(site_year, block, drop = TRUE),
    whole_plot_id = interaction(
      block_id, primary_method, drop = TRUE
    )
  )

validation_summary <- plot_data |>
  group_by(site_year) |>
  summarise(
    plots = n(),
    blocks = n_distinct(block_id),
    whole_plots = n_distinct(whole_plot_id),
    min_biomass = min(bean_biomass_kg_ha),
    max_biomass = max(bean_biomass_kg_ha),
    .groups = "drop"
  )

show_table(
  validation_summary,
  "Validated records and soybean biomass range (kg/ha)"
)
```

| site_year | plots | blocks | whole_plots | min_biomass | max_biomass |
|:----------|------:|-------:|------------:|------------:|------------:|
| V_2023    |    20 |      4 |           8 |      3997.4 |      6923.5 |
| O2E_2024  |    20 |      4 |           8 |      7159.3 |      9524.1 |
| O2W_2024  |    20 |      4 |           8 |      6446.0 |     10173.1 |

Validated records and soybean biomass range (kg/ha)

# Prepare subplot observations

Each subplot contributes one biomass observation, so no within-subplot
averaging is needed. Sum contrasts support Type III fixed-effect tests
in the model containing the management-system-by-site-year interaction.

``` r
biomass_data <- plot_data

contrasts(biomass_data$treatment) <- contr.sum(5)
contrasts(biomass_data$site_year) <- contr.sum(3)
```

# Describe the observations

The table summarizes observed means and standard deviations for each
management system within each site-year. Standard deviations describe
variation among the four subplot measurements and are not model-based
standard errors.

``` r
biomass_summary <- biomass_data |>
  group_by(site_year, treatment) |>
  summarise(
    n_plots = n(),
    mean_biomass = mean(bean_biomass_kg_ha),
    sd_biomass = sd(bean_biomass_kg_ha),
    .groups = "drop"
  )

show_table(
  biomass_summary,
  "Observed subplot soybean biomass: mean and SD (kg/ha)"
)
```

| site_year | treatment | n_plots | mean_biomass | sd_biomass |
|:----------|:----------|--------:|-------------:|-----------:|
| V_2023    | RNO       |       4 |       4731.6 |      509.9 |
| V_2023    | RIM       |       4 |       4552.4 |      353.6 |
| V_2023    | RIC       |       4 |       5167.6 |     1056.0 |
| V_2023    | TIM       |       4 |       5072.6 |     1072.4 |
| V_2023    | TIC       |       4 |       5637.2 |      883.0 |
| O2E_2024  | RNO       |       4 |       8364.9 |      741.7 |
| O2E_2024  | RIM       |       4 |       8061.0 |      445.9 |
| O2E_2024  | RIC       |       4 |       8444.7 |      939.2 |
| O2E_2024  | TIM       |       4 |       9086.5 |      429.0 |
| O2E_2024  | TIC       |       4 |       8878.4 |      120.2 |
| O2W_2024  | RNO       |       4 |       7976.3 |      755.6 |
| O2W_2024  | RIM       |       4 |       7593.6 |      426.1 |
| O2W_2024  | RIC       |       4 |       7468.9 |      964.6 |
| O2W_2024  | TIM       |       4 |       9261.5 |      754.3 |
| O2W_2024  | TIC       |       4 |       9544.9 |      682.1 |

Observed subplot soybean biomass: mean and SD (kg/ha)

The exploratory figure displays individual observations and arithmetic
means. It helps reveal site-year differences and differences in observed
spread before fitting the mixed model.

``` r
ggplot(biomass_data, aes(treatment, bean_biomass_kg_ha)) +
  geom_point(
    aes(color = treatment),
    position = position_jitter(width = 0.12, height = 0, seed = 123),
    size = 2.5,
    alpha = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4,
    color = "black"
  ) +
  facet_wrap(
    ~site_year,
    labeller = labeller(site_year = experiment_labels)
  ) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Management system",
    y = "Soybean crop biomass (kg/ha)",
    caption = paste(
      "Colored points are subplot observations;",
      "black diamonds are treatment means."
    )
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.caption = element_text(size = 10, hjust = 0)
  )
```

![](soybean_crop_biomass-1-_files/figure-gfm/observed-figure-1.png)<!-- -->

# Fit the primary model

The initial model includes block and termination whole-plot random
intercepts. Management system and site-year, including their
interaction, are fixed effects. The interaction is retained so that
comparisons can be estimated separately by site-year or averaged equally
across site-years, regardless of the interaction p-value. Site-year is
fixed, so inference describes these three observed experiments. Unique
block IDs nest blocks within site-years, and whole-plot IDs nest
termination whole plots within blocks. A whole-plot intercept models
correlation among its subplots. With one observation per subplot, a
separate subplot intercept cannot be distinguished from residual error.
The model assumes homogeneous residual variance.

``` r
biomass_model_full <- lmerTest::lmer(
  bean_biomass_kg_ha ~ treatment * site_year +
    (1 | block_id) + (1 | whole_plot_id),
  data = biomass_data,
  REML = TRUE
)
```

    ## boundary (singular) fit: see help('isSingular')

# Evaluate random-effects structure

The initial model estimated the additional block-level variance at zero
and reported a singular fit. We retain it for documentation and fit a
simpler model with the whole-plot intercept. The latter preserves
within-whole-plot correlation and is used for the subsequent analysis.
Singularity alone does not establish a treatment effect or its absence.

``` r
show_table(
  as.data.frame(lme4::VarCorr(biomass_model_full)) |>
    transmute(component = grp, variance = vcov, SD = sdcor),
  "Variance components for initial biomass model", digits = 1
)
```

| component     | variance |    SD |
|:--------------|---------:|------:|
| whole_plot_id | 215890.5 | 464.6 |
| block_id      |      0.0 |   0.0 |
| Residual      | 313585.5 | 560.0 |

Variance components for initial biomass model

``` r
cat("Initial model singular:", lme4::isSingular(biomass_model_full), "\n")
```

    ## Initial model singular: TRUE

``` r
biomass_model <- lmerTest::lmer(
  bean_biomass_kg_ha ~ treatment * site_year +
    (1 | whole_plot_id),
  data = biomass_data,
  REML = TRUE
)

cat(
  "Biomass model singular:",
  lme4::isSingular(biomass_model, tol = 1e-4),
  "\n"
)
```

    ## Biomass model singular: FALSE

# Sensitivity to the random-effects structure

The initial model includes a block intercept whose estimated variance
was zero. To check whether omitting it changes inference, compare the
treatment estimated marginal means and their SEs from both models,
averaging site-years equally. We also compare all pairwise treatment
contrasts. This is a sensitivity check, not a new model-selection test.
Subsequent results continue to use the simpler `biomass_model`.

``` r
full_emm <- emmeans(biomass_model_full, ~ treatment, weights = "equal")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
simple_emm <- emmeans(biomass_model, ~ treatment, weights = "equal")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
sensitivity_means <- as.data.frame(summary(full_emm)) |>
  select(treatment, full_mean = emmean, full_SE = SE, full_df = df) |>
  left_join(
    as.data.frame(summary(simple_emm)) |>
      select(treatment, simple_mean = emmean, simple_SE = SE, simple_df = df),
    by = "treatment"
  ) |>
  mutate(
    mean_difference = simple_mean - full_mean,
    SE_difference = simple_SE - full_SE
  )

show_table(
  sensitivity_means,
  "Treatment EMM sensitivity to the random-effects structure (kg/ha)",
  digits = 2
)
```

| treatment | full_mean | full_SE | full_df | simple_mean | simple_SE | simple_df | mean_difference | SE_difference |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| RNO | 7024.25 | 210.05 | 35.88 | 7024.25 | 210.05 | 35.88 | 0 | 0 |
| RIM | 6735.67 | 210.05 | 35.88 | 6735.67 | 210.05 | 35.88 | 0 | 0 |
| RIC | 7027.07 | 210.05 | 35.88 | 7027.07 | 210.05 | 35.88 | 0 | 0 |
| TIM | 7806.88 | 210.05 | 35.88 | 7806.88 | 210.05 | 35.88 | 0 | 0 |
| TIC | 8020.17 | 210.05 | 35.88 | 8020.17 | 210.05 | 35.88 | 0 | 0 |

Treatment EMM sensitivity to the random-effects structure (kg/ha)

``` r
full_pairs <- as.data.frame(summary(pairs(full_emm, adjust = "none"))) |>
  select(contrast, full_estimate = estimate, full_SE = SE,
         full_df = df, full_p = p.value)
simple_pairs <- as.data.frame(summary(pairs(simple_emm, adjust = "none"))) |>
  select(contrast, simple_estimate = estimate, simple_SE = SE,
         simple_df = df, simple_p = p.value)

sensitivity_pairs <- full_pairs |>
  left_join(simple_pairs, by = "contrast") |>
  mutate(
    estimate_difference = simple_estimate - full_estimate,
    SE_difference = simple_SE - full_SE
  )

show_table(
  sensitivity_pairs,
  "Pairwise treatment-contrast sensitivity (kg/ha; unadjusted p-values)",
  digits = 4
)
```

| contrast | full_estimate | full_SE | full_df | full_p | simple_estimate | simple_SE | simple_df | simple_p | estimate_difference | SE_difference |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| RNO - RIM | 288.5833 | 228.6137 | 27.3281 | 0.2175 | 288.5833 | 228.6137 | 27.3281 | 0.2175 | 0 | 0 |
| RNO - RIC | -2.8167 | 228.6137 | 27.3281 | 0.9903 | -2.8167 | 228.6137 | 27.3281 | 0.9903 | 0 | 0 |
| RNO - TIM | -782.6250 | 297.0623 | 35.8791 | 0.0124 | -782.6250 | 297.0623 | 35.8791 | 0.0124 | 0 | 0 |
| RNO - TIC | -995.9167 | 297.0623 | 35.8791 | 0.0019 | -995.9167 | 297.0623 | 35.8791 | 0.0019 | 0 | 0 |
| RIM - RIC | -291.4000 | 228.6137 | 27.3281 | 0.2132 | -291.4000 | 228.6137 | 27.3281 | 0.2132 | 0 | 0 |
| RIM - TIM | -1071.2083 | 297.0623 | 35.8791 | 0.0009 | -1071.2083 | 297.0623 | 35.8791 | 0.0009 | 0 | 0 |
| RIM - TIC | -1284.5000 | 297.0623 | 35.8791 | 0.0001 | -1284.5000 | 297.0623 | 35.8791 | 0.0001 | 0 | 0 |
| RIC - TIM | -779.8083 | 297.0623 | 35.8791 | 0.0126 | -779.8083 | 297.0623 | 35.8791 | 0.0126 | 0 | 0 |
| RIC - TIC | -993.1000 | 297.0623 | 35.8791 | 0.0019 | -993.1000 | 297.0623 | 35.8791 | 0.0019 | 0 | 0 |
| TIM - TIC | -213.2917 | 228.6137 | 27.3281 | 0.3590 | -213.2917 | 228.6137 | 27.3281 | 0.3590 | 0 | 0 |

Pairwise treatment-contrast sensitivity (kg/ha; unadjusted p-values)

``` r
# Compare both reporting procedures, using the same families in each model.
sensitivity_decisions <- tibble(
  procedure = c("Fisher's LSD (unadjusted)", "Tukey (ten comparisons)"),
  changed_decisions = vapply(c("none", "tukey"), function(adjustment) {
    full_p <- as.data.frame(summary(pairs(full_emm, adjust = adjustment)))$p.value
    simple_p <- as.data.frame(summary(pairs(simple_emm, adjust = adjustment)))$p.value
    sum((full_p < 0.05) != (simple_p < 0.05))
  }, integer(1))
)
show_table(sensitivity_decisions, "Changes in pairwise decisions at α = 0.05", digits = 0)
```

| procedure                 | changed_decisions |
|:--------------------------|------------------:|
| Fisher’s LSD (unadjusted) |                 0 |
| Tukey (ten comparisons)   |                 0 |

Changes in pairwise decisions at α = 0.05

The initial model’s estimated block variance was 0.0. Dropping that
intercept changes pooled treatment means by at most 9.09e-12 kg/ha and
pairwise contrast SEs by at most 1.19e-05 kg/ha. There are 0 changed
decisions across the two comparison procedures. These checks assess
numerical stability, not whether blocking was unnecessary in the design.
EMMs are compared using means, SEs, and degrees of freedom; p-values
test treatment differences rather than whether an individual biomass
mean differs from zero.

# Check model assumptions

Variance components quantify the remaining between-whole-plot and
residual variation **after accounting for the fixed effects**. They do
not measure the amount of biomass variation caused by termination
method. The plots assess residual spread and approximate normality; they
do not prove the assumptions.

``` r
variance_table <- as.data.frame(lme4::VarCorr(biomass_model)) |>
  transmute(
    component = grp,
    variance = vcov,
    SD = sdcor
  )

show_table(
  variance_table,
  "Variance components for soybean biomass (kg/ha)",
  digits = 1
)
```

| component     | variance |    SD |
|:--------------|---------:|------:|
| whole_plot_id | 215890.6 | 464.6 |
| Residual      | 313585.4 | 560.0 |

Variance components for soybean biomass (kg/ha)

``` r
cat(
  "Singular fit:",
  lme4::isSingular(biomass_model, tol = 1e-4),
  "\n"
)
```

    ## Singular fit: FALSE

``` r
diag_data <- biomass_data |>
  mutate(
    fitted = fitted(biomass_model),
    residual = residuals(biomass_model)
  )
```

``` r
ggplot(diag_data, aes(fitted, residual, color = site_year)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_point(size = 2.5) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = "Fitted soybean biomass (kg/ha)",
    y = "Residual (kg/ha)",
    color = "Site-year"
  )
```

![](soybean_crop_biomass-1-_files/figure-gfm/residual-diagnostics-1.png)<!-- -->

``` r
ggplot(diag_data, aes(sample = residual)) +
  stat_qq() +
  stat_qq_line() +
  labs(
    x = "Theoretical normal quantile",
    y = "Residual quantile"
  )
```

![](soybean_crop_biomass-1-_files/figure-gfm/residual-diagnostics-2.png)<!-- -->

``` r
ggplot(diag_data, aes(site_year, residual, color = treatment)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 123),
    size = 2.5
  ) +
  scale_color_manual(values = treatment_colors) +
  scale_x_discrete(labels = experiment_labels) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = NULL,
    y = "Residual (kg/ha)",
    color = "Management system"
  )
```

![](soybean_crop_biomass-1-_files/figure-gfm/residual-diagnostics-3.png)<!-- -->
\# Test management system, site-year, and their interaction

Type III tests use Satterthwaite degrees of freedom with sum contrasts.
The management-system test addresses mean differences averaged over the
three site-years; the interaction assesses evidence that those
differences vary across site-years.

``` r
biomass_anova <- anova(
  biomass_model,
  type = 3,
  ddf = "Satterthwaite"
)

show_table(
  as.data.frame(biomass_anova) |>
    tibble::rownames_to_column("term"),
  "Type III ANOVA for soybean biomass",
  digits = 3
)
```

| term                |   Sum Sq |    Mean Sq | NumDF |  DenDF | F value | Pr(\>F) |
|:--------------------|---------:|-----------:|------:|-------:|--------:|--------:|
| treatment           |  6208297 |  1552074.2 |     4 | 33.409 |   4.949 |   0.003 |
| site_year           | 56554818 | 28277409.1 |     2 | 17.484 |  90.174 |   0.000 |
| treatment:site_year |  3131650 |   391456.3 |     8 | 33.409 |   1.248 |   0.303 |

Type III ANOVA for soybean biomass

``` r
readr::write_csv(as.data.frame(biomass_anova) |> tibble::rownames_to_column("term"),
                 file.path(tab_dir, "primary_model_anova.csv"))
```

The management-system test gave F = 4.949 (p = 0.003), and the site-year
test gave F = 90.174 (p \< 0.001). The interaction gave F = 1.248 (p =
0.303). A nonsignificant interaction does not establish identical system
effects. The interaction remains in the model regardless of this test.

# Estimated marginal means

Model-based means average the three site-years equally. Because the
design is balanced, these may coincide closely with the observed
treatment means. Model-based SEs reflect the covariance structure and
should not be confused with the descriptive SEs calculated from raw
subplot observations.

``` r
biomass_emm <- simple_emm

biomass_emm_table <- as.data.frame(summary(biomass_emm)) |>
  select(treatment, emmean, SE, df, lower.CL, upper.CL)

show_table(
  biomass_emm_table,
  "Estimated marginal mean soybean biomass by management system (kg/ha)",
  digits = 1
)
```

| treatment | emmean |    SE |   df | lower.CL | upper.CL |
|:----------|-------:|------:|-----:|---------:|---------:|
| RNO       | 7024.3 | 210.1 | 35.9 |   6598.2 |   7450.3 |
| RIM       | 6735.7 | 210.1 | 35.9 |   6309.6 |   7161.7 |
| RIC       | 7027.1 | 210.1 | 35.9 |   6601.0 |   7453.1 |
| TIM       | 7806.9 | 210.1 | 35.9 |   7380.8 |   8232.9 |
| TIC       | 8020.2 | 210.1 | 35.9 |   7594.1 |   8446.2 |

Estimated marginal mean soybean biomass by management system (kg/ha)

``` r
readr::write_csv(biomass_emm_table,
                 file.path(tab_dir, "treatment_model_means.csv"))
```

The following estimates retain the interaction and describe each
site-year separately. SEs come from the primary common-variance mixed
model. No additional site-year-specific pairwise tests are introduced.

``` r
biomass_cell_emm <- emmeans(biomass_model, ~ treatment | site_year)
biomass_cell_table <- as.data.frame(summary(biomass_cell_emm)) |>
  mutate(site_year = unname(experiment_labels[as.character(site_year)]),
         `Estimated mean ± SE (kg/ha)` = sprintf("%.0f ± %.0f", emmean, SE)) |>
  select(site_year, treatment, `Estimated mean ± SE (kg/ha)`) |>
  pivot_wider(names_from = treatment, values_from = `Estimated mean ± SE (kg/ha)`)
show_table(biomass_cell_table, "Model-estimated treatment means ± SE within site-year (kg/ha)")
```

| site_year            | RNO        | RIM        | RIC        | TIM        | TIC        |
|:---------------------|:-----------|:-----------|:-----------|:-----------|:-----------|
| 2023 – Field V       | 4732 ± 364 | 4552 ± 364 | 5168 ± 364 | 5073 ± 364 | 5637 ± 364 |
| 2024 – Field O2 East | 8365 ± 364 | 8061 ± 364 | 8445 ± 364 | 9087 ± 364 | 8878 ± 364 |
| 2024 – Field O2 West | 7976 ± 364 | 7594 ± 364 | 7469 ± 364 | 9262 ± 364 | 9545 ± 364 |

Model-estimated treatment means ± SE within site-year (kg/ha)

``` r
readr::write_csv(as.data.frame(summary(biomass_cell_emm)),
                 file.path(tab_dir, "site_year_treatment_model_means.csv"))
```

# Compare management systems

Tukey-adjusted pairwise comparisons provide a family-wise
error-controlled reference for all ten comparisons among five management
systems. The split-plot model gives different standard errors for
within-termination and across-termination comparisons: whole-plot
variation cancels from a within-method difference but contributes to an
across-method difference. All comparisons describe implemented
management systems; they do not isolate termination method from
supplemental management in a fully crossed factorial.

``` r
biomass_pairs <- pairs(
  biomass_emm,
  adjust = "tukey"
)

biomass_pairs_table <- as.data.frame(summary(biomass_pairs)) |>
  select(contrast, estimate, SE, df, t.ratio, p.value)

show_table(
  biomass_pairs_table,
  "Pairwise comparisons among soybean management systems (kg/ha)",
  digits = 3
)
```

| contrast  |  estimate |      SE |     df | t.ratio | p.value |
|:----------|----------:|--------:|-------:|--------:|--------:|
| RNO - RIM |   288.583 | 228.614 | 27.328 |   1.262 |   0.716 |
| RNO - RIC |    -2.817 | 228.614 | 27.328 |  -0.012 |   1.000 |
| RNO - TIM |  -782.625 | 297.062 | 35.879 |  -2.635 |   0.085 |
| RNO - TIC |  -995.917 | 297.062 | 35.879 |  -3.353 |   0.015 |
| RIM - RIC |  -291.400 | 228.614 | 27.328 |  -1.275 |   0.708 |
| RIM - TIM | -1071.208 | 297.062 | 35.879 |  -3.606 |   0.008 |
| RIM - TIC | -1284.500 | 297.062 | 35.879 |  -4.324 |   0.001 |
| RIC - TIM |  -779.808 | 297.062 | 35.879 |  -2.625 |   0.087 |
| RIC - TIC |  -993.100 | 297.062 | 35.879 |  -3.343 |   0.016 |
| TIM - TIC |  -213.292 | 228.614 | 27.328 |  -0.933 |   0.882 |

Pairwise comparisons among soybean management systems (kg/ha)

``` r
readr::write_csv(biomass_pairs_table,
                 file.path(tab_dir, "treatment_comparisons_tukey.csv"))
```

The targeted contrasts compare supplemental management options within
the same termination method. Holm adjustment is applied to this
four-contrast family. These contrasts avoid conflating different
termination methods with supplemental weed management.

``` r
within_method <- contrast(
  biomass_emm,
  method = list(
    "TIM - TIC" = c(0, 0, 0, 1, -1),
    "RIM - RIC" = c(0, 1, -1, 0, 0),
    "RNO - RIM" = c(1, -1, 0, 0, 0),
    "RNO - RIC" = c(1, 0, -1, 0, 0)
  )
)

show_table(
  as.data.frame(summary(within_method, adjust = "holm")) |>
    select(contrast, estimate, SE, df, t.ratio, p.value),
  "Within-termination management comparisons (kg/ha)",
  digits = 3
)
```

| contrast  | estimate |      SE |     df | t.ratio | p.value |
|:----------|---------:|--------:|-------:|--------:|--------:|
| TIM - TIC | -213.292 | 228.614 | 27.328 |  -0.933 |   0.853 |
| RIM - RIC | -291.400 | 228.614 | 27.328 |  -1.275 |   0.853 |
| RNO - RIM |  288.583 | 228.614 | 27.328 |   1.262 |   0.853 |
| RNO - RIC |   -2.817 | 228.614 | 27.328 |  -0.012 |   0.990 |

Within-termination management comparisons (kg/ha)

``` r
readr::write_csv(as.data.frame(summary(within_method, adjust = "holm")) |>
                   select(contrast, estimate, SE, df, t.ratio, p.value),
                 file.path(tab_dir, "within_method_comparisons_holm.csv"))
```

# Fisher’s LSD and letter groups

Unadjusted Fisher’s LSD comparisons are reported for continuity with the
laboratory’s established reporting practice. The overall
management-system ANOVA is the protection step for the current LSD
interpretation; however, *protected* LSD still does not provide the same
family-wise error control as Tukey across ten pairwise comparisons. The
corresponding compact letter display (CLD) assigns `a` to the highest
estimated mean. The code computes unadjusted comparisons regardless of
the omnibus p-value; if that test ceases to be significant after a data
change, revisit the protected-LSD interpretation. Contrast-specific SEs
and degrees of freedom are used, rather than a single universal LSD
cutoff. Shared letters indicate that a difference was not detected at
the chosen threshold, **not** that treatments are equivalent.

``` r
# Fisher's LSD: unadjusted pairwise comparisons.
crop_biomass_lsd <- pairs(biomass_emm, adjust = "none")

show_table(
  as.data.frame(summary(crop_biomass_lsd)) |>
    select(contrast, estimate, SE, df, t.ratio, p.value),
  "Fisher's LSD comparisons among management systems (kg/ha)",
  digits = 3
)
```

| contrast  |  estimate |      SE |     df | t.ratio | p.value |
|:----------|----------:|--------:|-------:|--------:|--------:|
| RNO - RIM |   288.583 | 228.614 | 27.328 |   1.262 |   0.218 |
| RNO - RIC |    -2.817 | 228.614 | 27.328 |  -0.012 |   0.990 |
| RNO - TIM |  -782.625 | 297.062 | 35.879 |  -2.635 |   0.012 |
| RNO - TIC |  -995.917 | 297.062 | 35.879 |  -3.353 |   0.002 |
| RIM - RIC |  -291.400 | 228.614 | 27.328 |  -1.275 |   0.213 |
| RIM - TIM | -1071.208 | 297.062 | 35.879 |  -3.606 |   0.001 |
| RIM - TIC | -1284.500 | 297.062 | 35.879 |  -4.324 |   0.000 |
| RIC - TIM |  -779.808 | 297.062 | 35.879 |  -2.625 |   0.013 |
| RIC - TIC |  -993.100 | 297.062 | 35.879 |  -3.343 |   0.002 |
| TIM - TIC |  -213.292 | 228.614 | 27.328 |  -0.933 |   0.359 |

Fisher’s LSD comparisons among management systems (kg/ha)

``` r
# Compact letter display, ordered from highest to lowest biomass.
crop_biomass_cld <- multcomp::cld(
  biomass_emm,
  adjust = "none",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE
) |>
  as.data.frame() |>
  mutate(
    .group = gsub("\\s+", "", .group),
    treatment = factor(treatment, levels = treatment_levels)
  )

show_table(
  crop_biomass_cld |>
    select(treatment, emmean, SE, .group),
  "Fisher's LSD groupings for soybean crop biomass",
  digits = 1
)
```

|     | treatment | emmean |    SE | .group |
|:----|:----------|-------:|------:|:-------|
| 2   | TIC       | 8020.2 | 210.1 | a      |
| 1   | TIM       | 7806.9 | 210.1 | a      |
| 3   | RIC       | 7027.1 | 210.1 | b      |
| 4   | RNO       | 7024.3 | 210.1 | b      |
| 5   | RIM       | 6735.7 | 210.1 | b      |

Fisher’s LSD groupings for soybean crop biomass

``` r
readr::write_csv(as.data.frame(summary(crop_biomass_lsd)) |>
                   select(contrast, estimate, SE, df, t.ratio, p.value),
                 file.path(tab_dir, "treatment_comparisons_fisher_lsd.csv"))
```

# Model-estimated biomass figures

The first figure shows management-system estimated marginal means with
model-based SEs and model-derived Fisher’s LSD groupings. Its error bars
reflect precision under the fitted model, not the observed
within-treatment spread. The interaction remains in the model when means
are averaged over site-years.

``` r
crop_biomass_treatment_data <- crop_biomass_cld |>
  mutate(label_y = emmean + SE + 250)

crop_biomass_treatment_fig <- ggplot(
  crop_biomass_treatment_data,
  aes(x = treatment, y = emmean, fill = treatment)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6,
    color = "black"
  ) +
  scale_fill_manual(values = treatment_colors, guide = "none") +
  scale_x_discrete(labels = treatment_labels) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Management system",
    y = "Soybean crop biomass (kg/ha)",
    caption = paste("Model-estimated means, averaged equally across site-years; error bars: ±1 model-based SE.",
                    lsd_caption, sep = "\n")
  ) +
  theme_classic(base_size = 16) + figure_caption_theme

crop_biomass_treatment_fig
```

![](soybean_crop_biomass-1-_files/figure-gfm/crop-biomass-treatment-figure-1.png)<!-- -->

``` r
dir.create(here("figures", "soybean_crop_biomass"), showWarnings = FALSE, recursive = TRUE)
ggsave(
  here("figures", "soybean_crop_biomass", "treatment_model.png"),
  plot = crop_biomass_treatment_fig,
  width = 9, height = 6, dpi = 300, bg = "white"
)
```

# Alternative figures: observed means and descriptive SEs

These figures show arithmetic means and descriptive SEs calculated from
the observations. They reuse the **model-derived** Fisher’s LSD letters
rather than testing differences using raw SEs. The error bars and
letters thus represent different quantities. For a given management
system, the 12 observations come from distinct whole plots across three
site-years. Its descriptive SE pools site-year differences with
within-site-year variation and does not represent the precision of an
equally weighted, site-year-adjusted treatment mean.

``` r
crop_biomass_treatment_raw <- biomass_data |>
  group_by(treatment) |>
  summarise(
    n = n(),
    mean_biomass = mean(bean_biomass_kg_ha),
    sd_biomass = sd(bean_biomass_kg_ha),
    se_biomass = sd_biomass / sqrt(n),
    .groups = "drop"
  ) |>
  left_join(
    crop_biomass_cld |>
      select(treatment, .group),
    by = "treatment"
  ) |>
  mutate(label_y = mean_biomass + se_biomass + 250)

crop_biomass_treatment_raw_fig <- ggplot(
  crop_biomass_treatment_raw,
  aes(x = treatment, y = mean_biomass, fill = treatment)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6
  ) +
  scale_fill_manual(values = treatment_colors, guide = "none") +
  scale_x_discrete(labels = treatment_labels) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Management system",
    y = "Soybean crop biomass (kg/ha)",
    caption = paste("Observed arithmetic means; error bars: ±1 descriptive SE (SD/√12), pooling site-years.",
                    "Letters compare model-estimated means averaged equally across site-years.",
                    lsd_caption, sep = "\n")
  ) +
  theme_classic(base_size = 16) + figure_caption_theme

crop_biomass_treatment_raw_fig
```

![](soybean_crop_biomass-1-_files/figure-gfm/crop-biomass-treatment-raw-figure-1.png)<!-- -->

``` r
ggsave(
  here("figures", "soybean_crop_biomass", "treatment_raw.png"),
  plot = crop_biomass_treatment_raw_fig,
  width = 9, height = 6, dpi = 300, bg = "white"
)
```

# Site-year comparisons

Site-year estimated marginal means average equally across the five
management systems. Fisher’s LSD letters summarize model-based
comparisons among the three experiments; these should be understood as
comparisons among the observed site-years rather than an isolated causal
effect of calendar year.

``` r
crop_biomass_year_emm <- emmeans(
  biomass_model,
  ~ site_year,
  weights = "equal"
)
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
crop_biomass_year_lsd <- pairs(
  crop_biomass_year_emm,
  adjust = "none"
)

show_table(
  as.data.frame(summary(crop_biomass_year_lsd)) |>
    select(contrast, estimate, SE, df, t.ratio, p.value),
  "Fisher's LSD comparisons among site-years (kg/ha)",
  digits = 3
)
```

| contrast            |  estimate |      SE |     df | t.ratio | p.value |
|:--------------------|----------:|--------:|-------:|--------:|--------:|
| V_2023 - O2E_2024   | -3534.830 | 295.787 | 17.484 | -11.951 |   0.000 |
| V_2023 - O2W_2024   | -3336.745 | 295.787 | 17.484 | -11.281 |   0.000 |
| O2E_2024 - O2W_2024 |   198.085 | 295.787 | 17.484 |   0.670 |   0.512 |

Fisher’s LSD comparisons among site-years (kg/ha)

``` r
crop_biomass_year_data <- multcomp::cld(
  crop_biomass_year_emm,
  adjust = "none",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE
) |>
  as.data.frame() |>
  mutate(
    site_year = factor(site_year, levels = names(experiment_labels)),
    .group = gsub("\\s+", "", .group),
    label_y = emmean + SE + 250
  )

crop_biomass_year_fig <- ggplot(
  crop_biomass_year_data,
  aes(x = site_year, y = emmean, fill = site_year)
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6
  ) +
  scale_fill_manual(values = experiment_colors, guide = "none") +
  scale_x_discrete(labels = experiment_labels) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Site-year",
    y = "Soybean crop biomass (kg/ha)",
    caption = paste("Model-estimated means, averaged equally across five systems; error bars: ±1 model-based SE.",
                    lsd_caption, sep = "\n")
  ) +
  theme_classic(base_size = 16) + figure_caption_theme

crop_biomass_year_fig
```

![](soybean_crop_biomass-1-_files/figure-gfm/crop-biomass-site-year-figure-1.png)<!-- -->

``` r
ggsave(
  here("figures", "soybean_crop_biomass", "site_year_model.png"),
  plot = crop_biomass_year_fig,
  width = 9, height = 6, dpi = 300, bg = "white"
)

readr::write_csv(as.data.frame(summary(crop_biomass_year_emm)),
                 file.path(tab_dir, "site_year_model_means.csv"))
readr::write_csv(as.data.frame(summary(crop_biomass_year_lsd)) |>
                   select(contrast, estimate, SE, df, t.ratio, p.value),
                 file.path(tab_dir, "site_year_comparisons_fisher_lsd.csv"))
```

For comparison, the observed site-year figure uses the 20 subplot values
per experiment. Its descriptive SE ignores within-whole-plot
correlation, so it should not be interpreted as design-adjusted
precision for testing site-year differences.

``` r
crop_biomass_year_raw <- biomass_data |>
  group_by(site_year) |>
  summarise(
    n = n(),
    mean_biomass = mean(bean_biomass_kg_ha),
    sd_biomass = sd(bean_biomass_kg_ha),
    se_biomass = sd_biomass / sqrt(n),
    .groups = "drop"
  ) |>
  left_join(
    crop_biomass_year_data |>
      select(site_year, .group),
    by = "site_year"
  ) |>
  mutate(label_y = mean_biomass + se_biomass + 250)

crop_biomass_year_raw_fig <- ggplot(
  crop_biomass_year_raw,
  aes(x = site_year, y = mean_biomass, fill = site_year)
) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_text(
    aes(y = label_y, label = .group),
    size = 6
  ) +
  scale_fill_manual(values = experiment_colors, guide = "none") +
  scale_x_discrete(labels = experiment_labels) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Site-year",
    y = "Soybean crop biomass (kg/ha)",
    caption = paste("Observed arithmetic means; error bars: ±1 descriptive SE (SD/√20), ignoring whole-plot correlation.",
                    "Letters compare model-estimated site-year means averaged equally across five systems.",
                    lsd_caption, sep = "\n")
  ) +
  theme_classic(base_size = 16) + figure_caption_theme

crop_biomass_year_raw_fig
```

![](soybean_crop_biomass-1-_files/figure-gfm/crop-biomass-site-year-raw-figure-1.png)<!-- -->

``` r
ggsave(
  here("figures", "soybean_crop_biomass", "site_year_raw.png"),
  plot = crop_biomass_year_raw_fig,
  width = 9, height = 6, dpi = 300, bg = "white"
)
```

# Results summary

Across the three observed site-years, soybean crop biomass differed
among the five management systems (p = 0.003) and among site-years (p \<
0.001). There was no statistically significant treatment-by-site-year
interaction (p = 0.303), but the model retains it to estimate system
means within individual site-years.

Estimated means were RNO: 7,024; RIM: 6,736; RIC: 7,027; TIM: 7,807;
TIC: 8,020 kg/ha. The tilled systems generally had higher biomass than
the roller-crimped systems. Fisher’s LSD detected differences in 6 of
the ten pairwise comparisons, covering all six tilled–rolled pairs.
Tukey detected 4 differences; it did not detect TIM–RNO (p = 0.085) or
TIM–RIC (p = 0.087). The figures’ letters describe LSD rather than Tukey
results. None of the four within-termination supplemental-management
comparisons was significant after Holm adjustment (smallest p = 0.853).
These system comparisons do not isolate the causal contribution of
termination method because the five systems are not a fully crossed
factorial.

Biomass was lower in 2023 Field V than in either 2024 experiment; the
two 2024 site-years did not differ by LSD (p = 0.512). These are
comparisons among observed experiments, not isolated effects of calendar
year. The 2023 RNO subplot that inadvertently received cultivation
limits interpretation of RNO versus RIC: the primary comparison
estimates differences by randomized assignment, including this
implementation deviation. No subplot has been reclassified.

Differences in soybean establishment, competition from living cereal
rye, planting conditions, and residue effects are possible explanations,
not mechanisms demonstrated by this biomass analysis. Greater crop
biomass does not necessarily imply greater grain yield; interpret these
results alongside establishment, weed biomass, and grain-yield analyses.

# Appendix: sensitivity to residual variance

The primary model assumes one residual variance across site-years.
Conditional residual plots assess spread and normality but do not prove
these assumptions. This check retains the treatment-by-site-year
interaction and whole-plot random intercept while allowing residual SD
to differ among site-years. It does not change treatment assignments or
replace the primary model.

``` r
# Fit both versions with the same engine and fixed effects. Their REML
# likelihoods are comparable because only the covariance structure changes.
# Pass contrasts through nlme's argument, rather than factor attributes:
# this avoids contrast-attribute warnings when emmeans reconstructs its grid.
variance_data <- biomass_data
contrasts(variance_data$treatment) <- NULL
contrasts(variance_data$site_year) <- NULL
variance_equal <- nlme::lme(
  fixed = bean_biomass_kg_ha ~ treatment * site_year,
  random = ~ 1 | whole_plot_id,
  data = variance_data, method = "REML",
  contrasts = list(treatment = contr.sum(5), site_year = contr.sum(3))
)
variance_by_year <- nlme::lme(
  fixed = bean_biomass_kg_ha ~ treatment * site_year,
  random = ~ 1 | whole_plot_id,
  weights = nlme::varIdent(form = ~ 1 | site_year),
  data = variance_data, method = "REML",
  contrasts = list(treatment = contr.sum(5), site_year = contr.sum(3))
)
variance_comparison <- as.data.frame(anova(variance_equal, variance_by_year))
show_table(variance_comparison |>
             select(-call) |> tibble::rownames_to_column("model"),
           "Common versus site-year-specific residual variance (REML)", digits = 3)
```

| model            | Model |  df |     AIC |     BIC |   logLik | Test   | L.Ratio | p-value |
|:-----------------|------:|----:|--------:|--------:|---------:|:-------|--------:|--------:|
| variance_equal   |     1 |  17 | 790.528 | 821.241 | -378.264 |        |      NA |      NA |
| variance_by_year |     2 |  19 | 792.594 | 826.921 | -377.297 | 1 vs 2 |   1.933 |    0.38 |

Common versus site-year-specific residual variance (REML)

``` r
residual_sd_ratios <- coef(variance_by_year$modelStruct$varStruct,
                          unconstrained = FALSE, allCoef = TRUE)
show_table(
  tibble(site_year = unname(experiment_labels[names(residual_sd_ratios)]),
         residual_SD_ratio = unname(residual_sd_ratios)),
  "Residual SD ratios relative to 2023 Field V", digits = 3
)
```

| site_year            | residual_SD_ratio |
|:---------------------|------------------:|
| 2023 – Field V       |             1.000 |
| 2024 – Field O2 East |             0.682 |
| 2024 – Field O2 West |             0.701 |

Residual SD ratios relative to 2023 Field V

``` r
variance_equal_emm <- emmeans(variance_equal, ~ treatment, weights = "equal")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
variance_year_emm <- emmeans(variance_by_year, ~ treatment, weights = "equal")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
variance_sensitivity <- as.data.frame(summary(pairs(variance_equal_emm, adjust = "none"))) |>
  select(contrast, common_estimate = estimate, common_SE = SE,
         common_df = df, common_p = p.value) |>
  left_join(
    as.data.frame(summary(pairs(variance_year_emm, adjust = "none"))) |>
      select(contrast, unequal_estimate = estimate, unequal_SE = SE,
             unequal_df = df, unequal_p = p.value),
    by = "contrast"
  )
show_table(variance_sensitivity,
           "Pooled contrast sensitivity to residual variance (kg/ha; unadjusted p-values)",
           digits = 4)
```

| contrast | common_estimate | common_SE | common_df | common_p | unequal_estimate | unequal_SE | unequal_df | unequal_p |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| RNO - RIM | 288.5833 | 228.6143 | 24 | 0.2190 | 288.5833 | 229.8074 | 24 | 0.2213 |
| RNO - RIC | -2.8167 | 228.6143 | 24 | 0.9903 | -2.8167 | 229.8074 | 24 | 0.9903 |
| RNO - TIM | -782.6250 | 297.0618 | 24 | 0.0145 | -782.6250 | 296.0916 | 24 | 0.0142 |
| RNO - TIC | -995.9167 | 297.0618 | 24 | 0.0026 | -995.9167 | 296.0916 | 24 | 0.0026 |
| RIM - RIC | -291.4000 | 228.6143 | 24 | 0.2146 | -291.4000 | 229.8074 | 24 | 0.2170 |
| RIM - TIM | -1071.2083 | 297.0618 | 24 | 0.0014 | -1071.2083 | 296.0916 | 24 | 0.0014 |
| RIM - TIC | -1284.5000 | 297.0618 | 24 | 0.0002 | -1284.5000 | 296.0916 | 24 | 0.0002 |
| RIC - TIM | -779.8083 | 297.0618 | 24 | 0.0148 | -779.8083 | 296.0916 | 24 | 0.0146 |
| RIC - TIC | -993.1000 | 297.0618 | 24 | 0.0027 | -993.1000 | 296.0916 | 24 | 0.0026 |
| TIM - TIC | -213.2917 | 228.6143 | 24 | 0.3601 | -213.2917 | 229.8074 | 24 | 0.3626 |

Pooled contrast sensitivity to residual variance (kg/ha; unadjusted
p-values)

``` r
# Pooled SEs can remain stable while site-year-specific precision changes.
variance_cell_se <- bind_rows(
  as.data.frame(emmeans(variance_equal, ~ treatment | site_year)) |>
    mutate(model = "Common residual variance"),
  as.data.frame(emmeans(variance_by_year, ~ treatment | site_year)) |>
    mutate(model = "Residual variance by site-year")
) |>
  group_by(site_year, model) |>
  summarise(min_SE = min(SE), max_SE = max(SE), .groups = "drop") |>
  mutate(site_year = unname(experiment_labels[as.character(site_year)]))
show_table(variance_cell_se, "Range of treatment-mean SEs within each site-year (kg/ha)")
```

| site_year            | model                          | min_SE | max_SE |
|:---------------------|:-------------------------------|-------:|-------:|
| 2023 – Field V       | Common residual variance       |  363.8 |  363.8 |
| 2023 – Field V       | Residual variance by site-year |  416.9 |  416.9 |
| 2024 – Field O2 East | Common residual variance       |  363.8 |  363.8 |
| 2024 – Field O2 East | Residual variance by site-year |  329.8 |  329.8 |
| 2024 – Field O2 West | Common residual variance       |  363.8 |  363.8 |
| 2024 – Field O2 West | Residual variance by site-year |  334.6 |  334.6 |

Range of treatment-mean SEs within each site-year (kg/ha)

The residual-variance comparison gave p = 0.380; AIC was 790.53 for
common variance and 792.59 for unequal variance. There were 0 changed
pooled LSD decisions between these two sensitivity fits. Their
unadjusted p-values use nlme/emmeans containment degrees of freedom
rather than the primary model’s Satterthwaite approximation; compare
covariance assumptions within this pair of fits rather than attributing
every difference from the primary analysis to unequal variance.

These data provide limited evidence for different residual variances and
no compelling reason to replace the primary common-variance model. This
does not establish equal variances. Pooled comparisons can be stable
while uncertainty within individual site-years changes, as the SE table
shows. No additional site-year-specific pairwise comparisons are
performed.

# Reproducibility and remaining limitations

The processed-data checks preserve one originally labeled treatment per
block; they cannot independently establish the randomized identity of
the affected 2023 subplot without its field record. Any discrepancy with
that record should be flagged before analysis, not corrected by silently
changing treatment labels. The operations error does not change the
whole-plot grouping. Interpret RNO–RIC comparisons with the
implementation deviation in mind.

The analysis uses the processed RDS directly. Packages include
tidyverse, here, conflicted, WrensBookshelf, lmerTest/lme4, emmeans,
multcomp, multcompView, nlme, knitr, and rmarkdown; the session below
records their environment. Render from the repository with these
dependencies installed. Four final PNGs are exported to `figures/` using
the existing filenames and dimensions.

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Tahoe 26.7.1
    ## 
    ## Matrix products: default
    ## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] emmeans_1.11.2-8     lmerTest_3.1-3       lme4_1.1-37         
    ##  [4] Matrix_1.7-4         WrensBookshelf_0.1.0 conflicted_1.2.0    
    ##  [7] here_1.0.2           lubridate_1.9.4      forcats_1.0.1       
    ## [10] stringr_1.5.2        dplyr_1.1.4          purrr_1.1.0         
    ## [13] readr_2.1.5          tidyr_1.3.1          tibble_3.3.0        
    ## [16] ggplot2_4.0.0        tidyverse_2.0.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6        xfun_0.61           lattice_0.22-7     
    ##  [4] tzdb_0.5.0          numDeriv_2016.8-1.1 vctrs_0.6.5        
    ##  [7] tools_4.5.2         Rdpack_2.6.4        generics_0.1.4     
    ## [10] parallel_4.5.2      sandwich_3.1-1      pkgconfig_2.0.3    
    ## [13] RColorBrewer_1.1-3  S7_0.2.0            lifecycle_1.0.4    
    ## [16] compiler_4.5.2      farver_2.1.2        textshaping_1.0.4  
    ## [19] codetools_0.2-20    htmltools_0.5.8.1   yaml_2.3.10        
    ## [22] crayon_1.5.3        pillar_1.11.1       nloptr_2.2.1       
    ## [25] MASS_7.3-65         cachem_1.1.0        reformulas_0.4.4   
    ## [28] boot_1.3-32         multcomp_1.4-29     nlme_3.1-168       
    ## [31] tidyselect_1.2.1    digest_0.6.37       mvtnorm_1.3-3      
    ## [34] stringi_1.8.7       labeling_0.4.3      splines_4.5.2      
    ## [37] rprojroot_2.1.1     fastmap_1.2.0       grid_4.5.2         
    ## [40] cli_3.6.5           magrittr_2.0.4      survival_3.8-3     
    ## [43] TH.data_1.1-4       withr_3.0.2         scales_1.4.0       
    ## [46] bit64_4.6.0-1       timechange_0.3.0    estimability_1.5.1 
    ## [49] rmarkdown_2.30      bit_4.6.0           ragg_1.5.0         
    ## [52] zoo_1.8-14          hms_1.1.4           memoise_2.0.1      
    ## [55] coda_0.19-4.1       evaluate_1.0.5      knitr_1.50         
    ## [58] rbibutils_2.3       rlang_1.1.6         Rcpp_1.1.0         
    ## [61] xtable_1.8-4        glue_1.8.0          vroom_1.6.6        
    ## [64] rstudioapi_0.17.1   minqa_1.2.8         R6_2.6.1           
    ## [67] systemfonts_1.3.1   multcompView_0.1-10
