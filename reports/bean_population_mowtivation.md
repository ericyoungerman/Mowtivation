Soybean population
================

- [Packages](#packages)
- [Data import & prep](#data-import--prep)
- [Model testing](#model-testing)
  - [Exploratory Analysis: Soybean
    population](#exploratory-analysis-soybean-population)
  - [Post-hoc summary table](#post-hoc-summary-table)
  - [ANOVA-style summary tables for bean
    population](#anova-style-summary-tables-for-bean-population)
- [Figures](#figures)
  - [Pooled model](#pooled-model)
  - [Pooled Raw](#pooled-raw)

\#Setup

# Packages

``` r
# Packages
library(tidyverse) # includes dplyr, ggplot2, readr, tibble, etc.
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

# Handle conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

# Treatment level order (use everywhere)
mow_levels <- c(
"Rolled, no control",
"Rolled + mowing",
"Rolled + high-residue cult.",
"Tilled + mowing",
"Tilled + cultivation"
)

# One consistent color palette for all figures

fill_cols <- c(
"Rolled, no control" = "#0072B2", # blue
"Rolled + mowing" = "#009E73", # green
"Rolled + high-residue cult." = "#F0E442", # yellow
"Tilled + mowing" = "#D55E00", # reddish
"Tilled + cultivation" = "#CC79A7" # magenta
)

#x-axis label helpers

label_break_spaces <- function(x) {
  stringr::str_replace_all(x, " ", "\n")
}

label_break_plus <- function(x) {
  stringr::str_replace_all(x, " \\+ ", "\n+ ")
}


#Helper: tidy emmeans output regardless of CI column names
#(works directly on an emmeans object)

tidy_emm <- function(emm, ref_levels = NULL) {
emm_df <- as.data.frame(emm)

lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(emm_df))[1]
ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(emm_df))[1]

if (is.na(lcl_col) || is.na(ucl_col)) {
stop("Could not find CI columns in emmeans output.")
}

out <- emm_df |>
mutate(
ci_low = .data[[lcl_col]],
ci_high = .data[[ucl_col]]
)

if (!is.null(ref_levels) && "weed_trt" %in% names(out)) {
out <- out |>
mutate(weed_trt = factor(weed_trt, levels = ref_levels))
}

out
}
```

# Data import & prep

``` r
bean_population_clean <- read_excel(
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
      "RIM" = "Rolled + mowing",
      "RIC" = "Rolled + high-residue cult.",
      "TIM" = "Tilled + mowing",
      "TIC" = "Tilled + cultivation"
    ),
    weed_trt = factor(weed_trt, levels = mow_levels)
  ) |>
  # keep only rows with non-missing population
  filter(!is.na(bean_population)) |>
  # optional per-area metrics, analogous to emergence
  mutate(
    bean_population_two_meter = bean_population * 2,
    bean_population_hectare   = (bean_population / 0.762) * 10000,
    bean_population_acre      = bean_population_hectare / 2.471
  )

# 2023 Field V subset (for the boxplot later)
bean_population_field_v_2023 <- bean_population_clean |>
  filter(year == 2023, location == "field v")

# Quick check
kable(
  head(bean_population_clean),
  caption = "All site-years, cleaned (bean population)"
)
```

| id | location | year | weed_trt | block | plot | bean_emergence | bean_biomass | inrow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | seed_weight | site_year | bean_population_two_meter | bean_population_hectare | bean_population_acre |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|
| CU_B1_P101 | field v | 2023 | Tilled + mowing | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 17.1200 | 2023.field v | 69 | 452755.9 | 183227.8 |
| CU_B1_P102 | field v | 2023 | Tilled + cultivation | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 17.4750 | 2023.field v | 79 | 518372.7 | 209782.6 |
| CU_B1_P103 | field v | 2023 | Rolled + mowing | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 7.840 | 37.5 | 449.93 | 16.7525 | 2023.field v | 75 | 492126.0 | 199160.7 |
| CU_B1_P104 | field v | 2023 | Rolled, no control | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35.0 | 412.59 | 16.1450 | 2023.field v | 70 | 459317.6 | 185883.3 |
| CU_B1_P105 | field v | 2023 | Rolled + high-residue cult. | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39.0 | 473.79 | 17.0475 | 2023.field v | 78 | 511811.0 | 207127.1 |
| CU_B1_P201 | field v | 2023 | Rolled + high-residue cult. | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 17.1500 | 2023.field v | 67 | 439632.5 | 177916.9 |

All site-years, cleaned (bean population)

# Model testing

### Exploratory Analysis: Soybean population

``` r
# 1) Summary table: bean population by site-year × treatment
bean_population_clean |>
  group_by(site_year, weed_trt) |>
  summarise(
    n      = n(),
    mean   = mean(bean_population_hectare, na.rm = TRUE),
    median = median(bean_population_hectare, na.rm = TRUE),
    sd     = sd(bean_population_hectare, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(site_year, weed_trt) |>
  kable(
    digits  = 1,
    caption = "Bean population (plants/ha) by site-year × treatment"
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Bean population (plants/ha) by site-year × treatment
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

587270.3
</td>

<td style="text-align:right;">

587270.3
</td>

<td style="text-align:right;">

46243.2
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

597112.9
</td>

<td style="text-align:right;">

600393.7
</td>

<td style="text-align:right;">

94785.4
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

631561.7
</td>

<td style="text-align:right;">

613517.1
</td>

<td style="text-align:right;">

42986.1
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

562664.0
</td>

<td style="text-align:right;">

567585.3
</td>

<td style="text-align:right;">

35284.9
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

579068.2
</td>

<td style="text-align:right;">

600393.7
</td>

<td style="text-align:right;">

59509.0
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

598753.3
</td>

<td style="text-align:right;">

593832.0
</td>

<td style="text-align:right;">

28538.9
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

564304.5
</td>

<td style="text-align:right;">

577427.8
</td>

<td style="text-align:right;">

45460.7
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

575787.4
</td>

<td style="text-align:right;">

587270.3
</td>

<td style="text-align:right;">

39868.2
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

570866.1
</td>

<td style="text-align:right;">

597112.9
</td>

<td style="text-align:right;">

56952.0
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 west
</td>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

606955.4
</td>

<td style="text-align:right;">

633202.1
</td>

<td style="text-align:right;">

69131.6
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

446194.2
</td>

<td style="text-align:right;">

449475.1
</td>

<td style="text-align:right;">

14174.8
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

470800.5
</td>

<td style="text-align:right;">

475721.8
</td>

<td style="text-align:right;">

40226.5
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

460958.0
</td>

<td style="text-align:right;">

475721.8
</td>

<td style="text-align:right;">

77523.2
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

508530.2
</td>

<td style="text-align:right;">

505249.3
</td>

<td style="text-align:right;">

52902.0
</td>

</tr>

<tr>

<td style="text-align:left;">

2023.field v
</td>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

524934.4
</td>

<td style="text-align:right;">

531496.1
</td>

<td style="text-align:right;">

25129.3
</td>

</tr>

</tbody>

</table>

``` r
# 2) Faceted boxplot: all site-years
bean_population_clean |>
  ggplot(aes(x = weed_trt, y = bean_population_hectare, fill = weed_trt)) +
  geom_boxplot(
    outlier.shape = NA,
    width  = 0.55,
    color  = "black"
  ) +
  geom_jitter(
    width  = 0.12,
    height = 0,
    alpha  = 0.4,
    size   = 1.8,
    color  = "grey30"
  ) +
  facet_wrap(~ site_year, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_plus) +
  labs(
    x     = NULL,
    y     = "Bean population (plants/ha)",
    title = "Bean population by treatment across site-years"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text  = element_text(face = "bold")
  )
```

![](figs/analysis/bean_yield-unnamed-chunk-3-1.png)<!-- -->

``` r
# 3) Boxplot: 2023 – Field V only
bean_population_field_v_2023 |>
  ggplot(aes(x = weed_trt, y = bean_population_hectare, fill = weed_trt)) +
  geom_boxplot(
    outlier.shape = NA,
    width  = 0.55,
    color  = "black"
  ) +
  geom_jitter(
    width  = 0.12,
    height = 0,
    alpha  = 0.4,
    size   = 1.8,
    color  = "grey30"
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_plus) +
  labs(
    x     = NULL,
    y     = "Bean population (plants/ha)",
    title = "Bean population by treatment, 2023 – Field V"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10)
  )
```

![](figs/analysis/bean_yield-unnamed-chunk-3-2.png)<!-- --> \##
Selection

``` r
### Model testing / selection for bean population (plants/ha)

options(contrasts = c("contr.sum", "contr.poly"))

# Interaction model: weed_trt * site_year --------------------------------
pop_int <- lmer(
  bean_population_hectare ~ weed_trt * site_year + (1 | site_year:block),
  data = bean_population_clean
)

# Additive model: weed_trt + site_year -----------------------------------
pop_add <- lmer(
  bean_population_hectare ~ weed_trt + site_year + (1 | site_year:block),
  data = bean_population_clean
)

# Compare models (AIC + LRT) ---------------------------------------------
aic_population <- tibble(
  model = c(
    "Additive: weed_trt + site_year",
    "Interaction: weed_trt * site_year"
  ),
  AIC = c(AIC(pop_add), AIC(pop_int))
)

kable(
  aic_population,
  digits  = 1,
  caption = "Bean population (plants/ha): model comparison (additive vs interaction)"
)
```

| model                              |    AIC |
|:-----------------------------------|-------:|
| Additive: weed_trt + site_year     | 1345.4 |
| Interaction: weed_trt \* site_year | 1181.5 |

Bean population (plants/ha): model comparison (additive vs interaction)

``` r
anova(pop_add, pop_int)  # LRT: is the interaction worth keeping?
```

    ## Data: bean_population_clean
    ## Models:
    ## pop_add: bean_population_hectare ~ weed_trt + site_year + (1 | site_year:block)
    ## pop_int: bean_population_hectare ~ weed_trt * site_year + (1 | site_year:block)
    ##         npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
    ## pop_add    9 1487.5 1506.4 -734.75    1469.5                     
    ## pop_int   17 1491.5 1527.1 -728.74    1457.5 12.028  8     0.1499

``` r
# Choose simpler additive model unless interaction is clearly needed -----
# (this is the model used in all downstream emmeans/plots)
pop.lmer <- pop_add

# Diagnostics on chosen model --------------------------------------------
set.seed(123)
res_pop <- DHARMa::simulateResiduals(pop.lmer)
plot(res_pop)
```

![](figs/analysis/bean_yield-unnamed-chunk-4-1.png)<!-- -->

``` r
DHARMa::testDispersion(pop.lmer)
```

![](figs/analysis/bean_yield-unnamed-chunk-4-2.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.89294, p-value = 0.552
    ## alternative hypothesis: two.sided

``` r
car::Anova(pop.lmer, type = 3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: bean_population_hectare
    ##                 Chisq Df Pr(>Chisq)    
    ## (Intercept) 6376.4358  1  < 2.2e-16 ***
    ## weed_trt       2.0854  4       0.72    
    ## site_year     51.5806  2  6.301e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Post-hoc summary table

``` r
### Bean population (plants/ha) with Fisher's LSD CLDs

# Estimated marginal means for weed_trt
emm_pop <- emmeans(pop.lmer, ~ weed_trt)

# Tidy emmeans (adds ci_low / ci_high and enforces treatment order)
emm_pop_df <- tidy_emm(emm_pop, ref_levels = mow_levels) |>
  as_tibble()

# Compact letter display (Fisher's LSD, no adjustment; "a" = highest)
cld_pop <- cld(
  emm_pop,
  adjust  = "none",
  Letters = letters,
  sort    = TRUE,   # default, but explicit
  reversed = TRUE   # now 'a' goes to the highest group(s)
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# Join emmeans + CLDs and format for reporting
emm_pop_df |>
  left_join(cld_pop, by = "weed_trt") |>
  select(weed_trt, emmean, SE, ci_low, ci_high, .group) |>
  mutate(across(c(emmean, SE, ci_low, ci_high), ~ round(.x, 1))) |>
  kable(
    caption   = "Estimated bean population (plants/ha) with 95% CI and Fisher's LSD group letters",
    col.names = c("Treatment", "Mean", "SE", "Lower CI", "Upper CI", "Group")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Estimated bean population (plants/ha) with 95% CI and Fisher’s LSD group
letters
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

544072.6
</td>

<td style="text-align:right;">

15468.1
</td>

<td style="text-align:right;">

513042.8
</td>

<td style="text-align:right;">

575102.4
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

544072.6
</td>

<td style="text-align:right;">

15468.1
</td>

<td style="text-align:right;">

513042.8
</td>

<td style="text-align:right;">

575102.4
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

556102.4
</td>

<td style="text-align:right;">

15468.1
</td>

<td style="text-align:right;">

525072.6
</td>

<td style="text-align:right;">

587132.1
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

547353.5
</td>

<td style="text-align:right;">

15468.1
</td>

<td style="text-align:right;">

516323.7
</td>

<td style="text-align:right;">

578383.2
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

570319.3
</td>

<td style="text-align:right;">

15468.1
</td>

<td style="text-align:right;">

539289.6
</td>

<td style="text-align:right;">

601349.1
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

### ANOVA-style summary tables for bean population

``` r
## 1) P-value summary (Location, Treatment, Interaction) -----------------

# Type-III tests for additive model (weed_trt + site_year)
anova_pop <- Anova(pop.lmer, type = 3)

anova_pop_df <- anova_pop |>
  as.data.frame() |>
  tibble::rownames_to_column("Effect")

# LRT for interaction (additive vs interaction models)
anova_interaction_pop <- anova(pop_add, pop_int)

pvals_pop <- tibble(
  Effect = c("Location (site_year)", "Treatment (weed_trt)", "Location × Treatment"),
  p_raw  = c(
    anova_pop_df$`Pr(>Chisq)`[anova_pop_df$Effect == "site_year"],
    anova_pop_df$`Pr(>Chisq)`[anova_pop_df$Effect == "weed_trt"],
    anova_interaction_pop$`Pr(>Chisq)`[2]
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

pvals_pop |>
  kable(
    caption   = "Bean population (plants/ha): P-values for location, treatment, and interaction",
    col.names = c("Effect", "P-value")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Bean population (plants/ha): P-values for location, treatment, and
interaction
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

0.720
</td>

</tr>

<tr>

<td style="text-align:left;">

Location × Treatment
</td>

<td style="text-align:left;">

0.150
</td>

</tr>

</tbody>

</table>

``` r
## 2) Location block: site-year means (model + raw) ----------------------

# Model-based emmeans by site_year
emm_loc_pop <- emmeans(pop.lmer, ~ site_year)

emm_loc_pop_df <- tidy_emm(emm_loc_pop) |>
  as_tibble() |>
  mutate(
    site_year  = as.factor(site_year),
    model_mean = emmean
  ) |>
  select(site_year, model_mean)

# CLDs for site_year (a = highest)
cld_loc_pop <- cld(
  emm_loc_pop,
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
raw_loc_pop <- bean_population_clean |>
  group_by(site_year) |>
  summarise(
    raw_mean = mean(bean_population_hectare, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(site_year = as.factor(site_year))

loc_summary_pop <- emm_loc_pop_df |>
  left_join(cld_loc_pop, by = "site_year") |>
  left_join(raw_loc_pop, by = "site_year") |>
  mutate(
    model_mean = round(model_mean, 1),
    raw_mean   = round(raw_mean, 1),
    raw_CLD    = loc_CLD  # use same letters for model + raw
  ) |>
  arrange(site_year)

loc_summary_pop |>
  kable(
    caption   = "Bean population (plants/ha): location (site-year) means with CLDs",
    col.names = c("Site-year", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Bean population (plants/ha): location (site-year) means with CLDs
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

591535.4
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

591535.4
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

583333.3
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

583333.3
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

482283.5
</td>

<td style="text-align:left;">

b
</td>

<td style="text-align:right;">

482283.5
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
emm_pop <- emmeans(pop.lmer, ~ weed_trt)

emm_trt_pop_df <- tidy_emm(emm_pop, ref_levels = mow_levels) |>
  as_tibble() |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    model_mean = emmean
  ) |>
  select(weed_trt, model_mean)

# CLDs for treatment (a = highest)
cld_pop <- cld(
  emm_pop,
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
raw_trt_pop <- bean_population_clean |>
  group_by(weed_trt) |>
  summarise(
    raw_mean = mean(bean_population_hectare, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(weed_trt = factor(weed_trt, levels = mow_levels))

trt_summary_pop <- emm_trt_pop_df |>
  left_join(cld_pop, by = "weed_trt") |>
  left_join(raw_trt_pop, by = "weed_trt") |>
  mutate(
    model_mean = round(model_mean, 1),
    raw_mean   = round(raw_mean, 1),
    raw_CLD    = trt_CLD
  ) |>
  arrange(weed_trt)

trt_summary_pop |>
  kable(
    caption   = "Bean population (plants/ha): treatment means with CLDs",
    col.names = c("Treatment", "Model mean", "Model CLD", "Raw mean", "Raw CLD")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Bean population (plants/ha): treatment means with CLDs
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

544072.6
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

544072.6
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

544072.6
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

544072.6
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

556102.4
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

556102.4
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

547353.5
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

547353.5
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

570319.3
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

570319.3
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
emm_sy_pop <- emmeans(pop.lmer, ~ weed_trt | site_year)

emm_sy_pop_df <- tidy_emm(emm_sy_pop, ref_levels = mow_levels) |>
  as_tibble() |>
  mutate(
    weed_trt   = factor(weed_trt, levels = mow_levels),
    site_year  = as.factor(site_year),
    model_mean = emmean
  ) |>
  select(site_year, weed_trt, model_mean)

# CLDs within each site_year (a = highest within that site_year)
cld_sy_pop <- cld(
  emm_sy_pop,
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
raw_sy_pop <- bean_population_clean |>
  group_by(site_year, weed_trt) |>
  summarise(
    raw_mean = mean(bean_population_hectare, na.rm = TRUE),
    .groups  = "drop"
  ) |>
  mutate(
    site_year = as.factor(site_year),
    weed_trt  = factor(weed_trt, levels = mow_levels)
  )

int_summary_pop <- emm_sy_pop_df |>
  left_join(cld_sy_pop, by = c("site_year", "weed_trt")) |>
  left_join(raw_sy_pop, by = c("site_year", "weed_trt")) |>
  mutate(
    model_mean = round(model_mean, 1),
    raw_mean   = round(raw_mean, 1),
    raw_CLD    = int_CLD
  ) |>
  arrange(site_year, weed_trt)

int_summary_pop |>
  kable(
    caption   = "Bean population (plants/ha): site-year × treatment means with CLDs",
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

Bean population (plants/ha): site-year × treatment means with CLDs
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

583224.0
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

587270.3
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled + mowing
</td>

<td style="text-align:right;">

583224.0
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

597112.9
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

595253.7
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

631561.7
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled + mowing
</td>

<td style="text-align:right;">

586504.8
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

562664.0
</td>

<td style="text-align:left;">

a
</td>

</tr>

<tr>

<td style="text-align:left;">

2024.field O2 east
</td>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

609470.7
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

579068.2
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

575021.9
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

598753.3
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

Rolled + mowing
</td>

<td style="text-align:right;">

575021.9
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

564304.5
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

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

587051.6
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

575787.4
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

Tilled + mowing
</td>

<td style="text-align:right;">

578302.7
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

570866.1
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

Tilled + cultivation
</td>

<td style="text-align:right;">

601268.6
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

606955.4
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

473972.0
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

446194.2
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

Rolled + mowing
</td>

<td style="text-align:right;">

473972.0
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

470800.5
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

Rolled + high-residue cult.
</td>

<td style="text-align:right;">

486001.7
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

460958.0
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

Tilled + mowing
</td>

<td style="text-align:right;">

477252.8
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

508530.2
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

Tilled + cultivation
</td>

<td style="text-align:right;">

500218.7
</td>

<td style="text-align:left;">

a
</td>

<td style="text-align:right;">

524934.4
</td>

<td style="text-align:left;">

a
</td>

</tr>

</tbody>

</table>

# Figures

## Pooled model

``` r
# Figure: Bean population by weed management treatment (pooled across site-years)

# Estimated marginal means for treatments
emm_pop <- emmeans(pop.lmer, ~ weed_trt)

# Tidy emmeans for plotting
emm_pop_df <- tidy_emm(emm_pop, ref_levels = mow_levels) |>
  as_tibble() |>
  mutate(
    response = emmean,
    ymin     = pmax(response - SE, 0),
    ymax     = response + SE
  )

# CLDs for treatments (Fisher's LSD, "a" = highest)
cld_pop <- cld(
  emm_pop,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# Merge means and CLDs
plot_df_pop <- emm_pop_df |>
  left_join(cld_pop, by = "weed_trt")

# Plot
ggplot(plot_df_pop, aes(x = weed_trt, y = response, fill = weed_trt)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust    = 0,
    fontface = "bold",
    size     = 6
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_plus) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x        = NULL,
    y        = "Bean population (plants/ha)",
    title    = "Bean population by weed management treatment",
    caption  = "Model-based means ± SE; letters = Fisher-style CLD for treatment main effect."
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title   = element_text(face = "bold")
  )
```

![](figs/analysis/bean_yield-unnamed-chunk-7-1.png)<!-- -->

``` r
# Save figure
ggsave(
  filename = here("figs", "analysis", "fig_bean_population_mowing_pooled.png"),
  width    = 7.5,
  height   = 5.5,
  dpi      = 300
)
```

## Pooled Raw

``` r
# Figure: Bean population by weed management treatment (raw means ± SE, model CLDs)

# 1) Raw means and SE by treatment --------------------------------------
raw_pop_summary <- bean_population_clean |>
  group_by(weed_trt) |>
  summarise(
    n    = n(),
    mean = mean(bean_population_hectare, na.rm = TRUE),
    sd   = sd(bean_population_hectare, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    ymin     = pmax(mean - se, 0),
    ymax     = mean + se
  )

# 2) Model-based CLDs for treatment main effect -------------------------
emm_pop <- emmeans(pop.lmer, ~ weed_trt)

cld_pop <- cld(
  emm_pop,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE   # "a" = highest mean
) |>
  as_tibble() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = str_trim(.group)
  ) |>
  select(weed_trt, .group)

# 3) Join raw means and model CLDs --------------------------------------
plot_df_pop_raw <- raw_pop_summary |>
  left_join(cld_pop, by = "weed_trt")

# 4) Plot ---------------------------------------------------------------
ggplot(plot_df_pop_raw, aes(x = weed_trt, y = mean, fill = weed_trt)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust    = 0,
    fontface = "bold",
    size     = 6
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = label_break_plus) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x        = NULL,
    y        = "Bean population (plants/ha)",
    title    = "Bean population by weed management treatment",
    caption  = "Raw means ± SE; letters = model-based Fisher-style CLD for treatment main effect."
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title   = element_text(face = "bold")
  )
```

![](figs/analysis/bean_yield-unnamed-chunk-8-1.png)<!-- -->

``` r
# Save figure
ggsave(
  filename = here("figs", "analysis", "fig_bean_population_mowing_pooled_raw.png"),
  width    = 7.5,
  height   = 5.5,
  dpi      = 300
)
```
