Soybean emergence
================

- [Data imort and prep](#data-imort-and-prep)

\#Setup

\#Packages

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
# Handle conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::recode)

# Helper: tidy emmeans regardless of CI column names
tidy_emm <- function(emm_df, ref_levels = NULL) {
  lcl_col <- intersect(c("lower.CL", "asymp.LCL"), names(emm_df))[1]
  ucl_col <- intersect(c("upper.CL", "asymp.UCL"), names(emm_df))[1]
  if (is.na(lcl_col) || is.na(ucl_col)) {
    stop("Could not find CI columns in emmeans output.")
  }
  out <- emm_df |>
    mutate(
      ci_low  = .data[[lcl_col]],
      ci_high = .data[[ucl_col]]
    )
  if (!is.null(ref_levels)) {
    out <- out |>
      mutate(weed_trt = factor(weed_trt, levels = ref_levels))
  }
  out
}
```

# Data imort and prep

``` r
# 1) Read + clean master
bean_emergence_clean <- read_excel(here("data", "raw", "All Treatments", "combined_raw.xlsx")) |>
  clean_names() |>
  rename(weed_trt = treatment) |>
  mutate(
    weed_trt = recode(
      weed_trt,
      "RNO" = "Rolled, no control",
      "RIM" = "Rolled + mowing",
      "RIC" = "Rolled + high-residue cult.",
      "TIM" = "Tilled + mowing",
      "TIC" = "Tilled + cultivation"
    ),
    site_year = factor(interaction(year, location, drop = TRUE)),
    block = factor(block),
    bean_emergence_two_meter = bean_emergence * 2,
    bean_emergence_acre = ((bean_emergence / 0.762) * 10000) / 2.471,
    bean_emergence_hectare = (bean_emergence / 0.762) * 10000
  ) |>
  filter(!is.na(bean_yield))

# 2) Lock in agronomic order
mow_levels <- c(
  "Rolled, no control",
  "Rolled + mowing",
  "Rolled + high-residue cult.",
  "Tilled + mowing",
  "Tilled + cultivation"
)

bean_emergence_clean <- bean_emergence_clean |>
  mutate(weed_trt = factor(weed_trt, levels = mow_levels))

# 3) Make the 2023 Field V subset
bean_emergence_field_v_2023 <- bean_emergence_clean |>
  filter(year == 2023, location == "field v") |>
  mutate(weed_trt = factor(weed_trt, levels = mow_levels))

# Quick check
kable(head(bean_emergence_clean), caption = "All site-years, cleaned (bean emergence)")
```

| id | location | year | weed_trt | block | plot | bean_emergence | bean_biomass | inrow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | seed_weight | site_year | bean_emergence_two_meter | bean_emergence_acre | bean_emergence_hectare |
|:---|:---|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|---:|
| CU_B1_P101 | field v | 2023 | Tilled + mowing | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 17.1200 | 2023.field v | 93 | 246959.2 | 610236.2 |
| CU_B1_P102 | field v | 2023 | Tilled + cultivation | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 17.4750 | 2023.field v | 85 | 225715.4 | 557742.8 |
| CU_B1_P103 | field v | 2023 | Rolled + mowing | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 7.840 | 37.5 | 449.93 | 16.7525 | 2023.field v | 73 | 193849.7 | 479002.6 |
| CU_B1_P104 | field v | 2023 | Rolled, no control | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35.0 | 412.59 | 16.1450 | 2023.field v | 82 | 217749.0 | 538057.7 |
| CU_B1_P105 | field v | 2023 | Rolled + high-residue cult. | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39.0 | 473.79 | 17.0475 | 2023.field v | 82 | 217749.0 | 538057.7 |
| CU_B1_P201 | field v | 2023 | Rolled + high-residue cult. | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 17.1500 | 2023.field v | 73 | 193849.7 | 479002.6 |

All site-years, cleaned (bean emergence)

\#Model testing \### Exploratory Analysis: Bean Emergence

``` r
# 1) Summary table: bean emergence by site-year × treatment
bean_emergence_clean |>
  group_by(site_year, weed_trt) |>
  summarise(
    n      = n(),
    mean   = mean(bean_emergence_hectare, na.rm = TRUE),
    median = median(bean_emergence_hectare, na.rm = TRUE),
    sd     = sd(bean_emergence_hectare, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(site_year, weed_trt) |>
  kable(
    digits  = 1,
    caption = "Bean emergence (plants/ha) by site-year × treatment"
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Bean emergence (plants/ha) by site-year × treatment
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

561023.6
</td>

<td style="text-align:right;">

577427.8
</td>

<td style="text-align:right;">

47165.0
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

585629.9
</td>

<td style="text-align:right;">

564304.5
</td>

<td style="text-align:right;">

47127.0
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

554461.9
</td>

<td style="text-align:right;">

551181.1
</td>

<td style="text-align:right;">

44341.9
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

588910.8
</td>

<td style="text-align:right;">

600393.7
</td>

<td style="text-align:right;">

81845.8
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

639763.8
</td>

<td style="text-align:right;">

636482.9
</td>

<td style="text-align:right;">

29099.2
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

554461.9
</td>

<td style="text-align:right;">

561023.6
</td>

<td style="text-align:right;">

20401.1
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

605315.0
</td>

<td style="text-align:right;">

583989.5
</td>

<td style="text-align:right;">

71353.6
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

600393.7
</td>

<td style="text-align:right;">

587270.3
</td>

<td style="text-align:right;">

79556.1
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

626640.4
</td>

<td style="text-align:right;">

623359.6
</td>

<td style="text-align:right;">

21762.6
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

593832.0
</td>

<td style="text-align:right;">

603674.5
</td>

<td style="text-align:right;">

35739.6
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

506889.8
</td>

<td style="text-align:right;">

508530.2
</td>

<td style="text-align:right;">

26988.1
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

498687.7
</td>

<td style="text-align:right;">

498687.7
</td>

<td style="text-align:right;">

26787.9
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

516732.3
</td>

<td style="text-align:right;">

508530.2
</td>

<td style="text-align:right;">

45578.9
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

590551.2
</td>

<td style="text-align:right;">

593832.0
</td>

<td style="text-align:right;">

20749.9
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

572506.6
</td>

<td style="text-align:right;">

577427.8
</td>

<td style="text-align:right;">

29527.6
</td>

</tr>

</tbody>

</table>

``` r
# 2) Faceted boxplot: all site-years
ggplot(bean_emergence_clean, aes(x = weed_trt, y = bean_emergence_hectare)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "black", fill = "#69b3a2") +
  geom_jitter(width = 0.12, height = 0, alpha = 0.4, size = 1.8, color = "#333333") +
  facet_wrap(~ site_year, nrow = 1, scales = "free_y") +
  scale_x_discrete(labels = \(x) stringr::str_replace_all(x, " ", "\n")) +
  labs(
    x     = NULL,
    y     = "Bean emergence (plants/ha)",
    title = "Bean emergence by treatment across site-years"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text  = element_text(face = "bold")
  )
```

![](figs/analysis/bean_emergence-unnamed-chunk-3-1.png)<!-- -->

``` r
# 3) Boxplot: 2023 – Field V only
bean_emergence_clean |>
  filter(year == 2023, location == "field v") |>
  ggplot(aes(x = weed_trt, y = bean_emergence_hectare)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "black", fill = "#69b3a2") +
  geom_jitter(width = 0.12, height = 0, alpha = 0.4, size = 1.8, color = "#333333") +
  scale_x_discrete(labels = \(x) stringr::str_replace_all(x, " ", "\n")) +
  labs(
    x     = NULL,
    y     = "Bean emergence (plants/ha)",
    title = "Bean emergence by treatment, 2023 – Field V"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(size = 10))
```

![](figs/analysis/bean_emergence-unnamed-chunk-3-2.png)<!-- -->
\##Selection

``` r
### Model testing / selection for bean emergence (plants/ha)

options(contrasts = c("contr.sum", "contr.poly"))

# Interaction model: weed_trt * site_year --------------------------------
emergence_int <- lmer(
  bean_emergence_hectare ~ weed_trt * site_year + (1 | site_year:block),
  data = bean_emergence_clean
)

# Additive model: weed_trt + site_year -----------------------------------
emergence_add <- lmer(
  bean_emergence_hectare ~ weed_trt + site_year + (1 | site_year:block),
  data = bean_emergence_clean
)

# Compare models (AIC + LRT) ---------------------------------------------
aic_emergence <- tibble(
  model = c("Additive: weed_trt + site_year",
            "Interaction: weed_trt * site_year"),
  AIC   = c(AIC(emergence_add), AIC(emergence_int))
)

kable(aic_emergence, digits = 1,
      caption = "Bean emergence (plants/ha): model comparison (additive vs interaction)")
```

| model                              |    AIC |
|:-----------------------------------|-------:|
| Additive: weed_trt + site_year     | 1331.4 |
| Interaction: weed_trt \* site_year | 1170.1 |

Bean emergence (plants/ha): model comparison (additive vs interaction)

``` r
anova(emergence_add, emergence_int)  # LRT: is the interaction worth keeping?
```

    ## Data: bean_emergence_clean
    ## Models:
    ## emergence_add: bean_emergence_hectare ~ weed_trt + site_year + (1 | site_year:block)
    ## emergence_int: bean_emergence_hectare ~ weed_trt * site_year + (1 | site_year:block)
    ##               npar    AIC    BIC  logLik -2*log(L) Chisq Df Pr(>Chisq)
    ## emergence_add    9 1471.8 1490.6 -726.90    1453.8                    
    ## emergence_int   17 1476.3 1511.9 -721.15    1442.3  11.5  8      0.175

``` r
# Choose simpler additive model unless interaction is clearly needed -----
emergence.lmer <- emergence_add

# Diagnostics on chosen model --------------------------------------------
res_emergence <- DHARMa::simulateResiduals(emergence.lmer)
plot(res_emergence)
```

![](figs/analysis/bean_emergence-unnamed-chunk-4-1.png)<!-- -->

``` r
DHARMa::testDispersion(emergence.lmer)
```

![](figs/analysis/bean_emergence-unnamed-chunk-4-2.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.88979, p-value = 0.528
    ## alternative hypothesis: two.sided

``` r
car::Anova(emergence.lmer, type = 3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: bean_emergence_hectare
    ##                Chisq Df Pr(>Chisq)    
    ## (Intercept) 7631.351  1  < 2.2e-16 ***
    ## weed_trt      17.262  4  0.0017194 ** 
    ## site_year     15.443  2  0.0004433 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

\###Post-hoc summary table

``` r
###  for bean emergence (plants/ha) with Fisher's LSD CLDs

# Estimated marginal means for weed_trt
emm_emergence <- emmeans(emergence.lmer, ~ weed_trt)

# Compact letter display (Fisher's LSD, no adjustment)
cld_emergence <- cld(emm_emergence, adjust = "none", Letters = letters)

# Tidy and format for reporting
cld_emergence |>
  as_tibble() |>
  select(weed_trt, emmean, SE, lower.CL, upper.CL, .group) |>
  mutate(across(c(emmean, SE, lower.CL, upper.CL), round, 1)) |>
  kable(
    caption = "Estimated means (plants/ha) with 95% CI and Fisher's LSD group letters",
    col.names = c("Treatment", "Mean", "SE", "Lower CI", "Upper CI", "Group")
  ) |>
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
```

<table class="table table-striped table-hover" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Estimated means (plants/ha) with 95% CI and Fisher’s LSD group letters
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

540791.8
</td>

<td style="text-align:right;">

13611.1
</td>

<td style="text-align:right;">

513473.6
</td>

<td style="text-align:right;">

568109.9
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

557196.0
</td>

<td style="text-align:right;">

13611.1
</td>

<td style="text-align:right;">

529877.8
</td>

<td style="text-align:right;">

584514.1
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

563210.8
</td>

<td style="text-align:right;">

13611.1
</td>

<td style="text-align:right;">

535892.7
</td>

<td style="text-align:right;">

590529.0
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

602034.1
</td>

<td style="text-align:right;">

13611.1
</td>

<td style="text-align:right;">

574716.0
</td>

<td style="text-align:right;">

629352.3
</td>

<td style="text-align:left;">

b
</td>

</tr>

<tr>

<td style="text-align:left;">

Tilled + cultivation
</td>

<td style="text-align:right;">

602034.1
</td>

<td style="text-align:right;">

13611.1
</td>

<td style="text-align:right;">

574716.0
</td>

<td style="text-align:right;">

629352.3
</td>

<td style="text-align:left;">

b
</td>

</tr>

</tbody>

</table>

\#Figures

``` r
#Figure: Bean emergence by weed management treatment (pooled across site-years)

# Estimated marginal means for treatments
emm_emergence <- emmeans(emergence.lmer, ~ weed_trt)

# Convert to data frame for plotting
emm_emergence_df <- emm_emergence |>
  as.data.frame() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    response = emmean,
    ymin     = pmax(response - SE, 0),
    ymax     = response + SE
  )

# CLDs for treatments (Fisher's LSD)
cld_emergence <- cld(
  emm_emergence,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as.data.frame() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = trimws(.group)
  )

# Merge means and CLDs
plot_df_emergence <- left_join(
  emm_emergence_df,
  cld_emergence |> select(weed_trt, .group),
  by = "weed_trt"
)

# Define colors (optional)
fill_cols <- viridis::viridis(length(mow_levels), option = "D", end = 0.9, begin = 0.1)

# Plot
ggplot(plot_df_emergence, aes(x = weed_trt, y = response, fill = weed_trt)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust = 0, fontface = "bold", size = 6
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = \(x) stringr::str_replace_all(x, " \\+ ", "\n+ ")) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    x = NULL,
    y = "Bean emergence (plants/ha)",
    title = "Bean emergence by weed management treatment (pooled across site-years)",
    caption = "Model-based means ± SE; letters = Fisher-style CLD for treatment main effect."
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    plot.title   = element_text(face = "bold")
  )
```

![](figs/analysis/bean_emergence-unnamed-chunk-6-1.png)<!-- -->

``` r
# Save figure
ggsave(
  filename = here::here("figs", "analysis", "fig_bean_emergence_mowing_pooled.png"),
  width    = 7.5,
  height   = 5.5,
  dpi      = 300
)
```

``` r
### Figure: Bean emergence by weed management treatment, faceted by site-year

# Estimated marginal means for weed_trt × site_year
emm_emergence_sy <- emmeans(emergence.lmer, ~ weed_trt | site_year)

# Convert to data frame for plotting
emm_sy_df <- emm_emergence_sy |>
  as.data.frame() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    response = emmean,
    ymin     = pmax(response - SE, 0),
    ymax     = response + SE
  )

# CLDs for each site-year
cld_sy <- cld(
  emm_emergence_sy,
  adjust   = "none",
  Letters  = letters,
  sort     = TRUE,
  reversed = TRUE
) |>
  as.data.frame() |>
  mutate(
    weed_trt = factor(weed_trt, levels = mow_levels),
    .group   = trimws(.group)
  )

# Merge means and CLDs
plot_df_sy <- left_join(
  emm_sy_df,
  cld_sy |> select(weed_trt, site_year, .group),
  by = c("weed_trt", "site_year")
)

# Define colors
fill_cols <- viridis::viridis(length(mow_levels), option = "D", end = 0.9, begin = 0.1)

# Plot faceted by site-year
ggplot(plot_df_sy, aes(x = weed_trt, y = response, fill = weed_trt)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.14) +
  geom_text(
    aes(y = ymax * 1.08, label = .group),
    vjust = 0, fontface = "bold", size = 5
  ) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = \(x) stringr::str_replace_all(x, " \\+ ", "\n+ ")) +
  scale_y_continuous(labels = scales::label_comma()) +  
  labs(
    x = NULL,
    y = "Bean emergence (plants/ha)",
    title = "Bean emergence by treatment within each site-year",
    caption = "Model-based means ± SE; letters = Fisher-style CLD for treatment effect within site-year."
  ) +
  facet_wrap(~ site_year, scales = "free_y") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x  = element_text(lineheight = 0.95, margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    strip.text   = element_text(face = "bold", size = 14),
    plot.title   = element_text(face = "bold")
  )
```

![](figs/analysis/bean_emergence-unnamed-chunk-7-1.png)<!-- -->

``` r
# Save figure
ggsave(
  filename = here::here("figs", "analysis", "fig_bean_emergence_by_site_year.png"),
  width    = 10,
  height   = 6.5,
  dpi      = 300
)
```
