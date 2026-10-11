Mowtivation Data Preparation
================

- [Setup](#setup)
- [Packages](#packages)
- [File paths](#file-paths)
- [Import raw data](#import-raw-data)
- [Initial data checks](#initial-data-checks)
  - [Site-year structure](#site-year-structure)
  - [Identifier checks](#identifier-checks)
  - [Missing values](#missing-values)
- [Prepare dataset](#prepare-dataset)
- [Conversions](#conversions)
- [Save procecessed data](#save-procecessed-data)

This document imports the raw Mowtivation plot-level dataset, checks its
structure, applies documented data preparation steps, and saves the
processed dataset used by the outcome-specific analyses.

# Setup

# Packages

``` r
library(tidyverse)
library(readxl)
library(janitor)
library(here)
```

# File paths

``` r
raw_file <- here(
  "data", "raw", "All Treatments", "combined_raw.xlsx"
)

processed_dir <- here("data", "processed")

stopifnot(file.exists(raw_file))
```

# Import raw data

``` r
raw_data <- read_excel(raw_file) |>
  clean_names()

glimpse(raw_data)
```

    ## Rows: 60
    ## Columns: 18
    ## $ id                     <chr> "V_2023_B1_P101", "V_2023_B1_P102", "V_2023_B1_…
    ## $ location               <chr> "V", "V", "V", "V", "V", "V", "V", "V", "V", "V…
    ## $ year                   <dbl> 2023, 2023, 2023, 2023, 2023, 2023, 2023, 2023,…
    ## $ block                  <dbl> 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4,…
    ## $ plot                   <dbl> 101, 102, 103, 104, 105, 201, 202, 203, 204, 20…
    ## $ treatment              <chr> "TIM", "TIC", "RIM", "RNO", "RIC", "RIC", "RNO"…
    ## $ cereal_rye_variety     <chr> "Aroostook", "Aroostook", "Aroostook", "Aroosto…
    ## $ bean_emergence         <dbl> 46.5, 42.5, 36.5, 41.0, 41.0, 36.5, 39.0, 40.0,…
    ## $ bean_biomass           <dbl> 223.740, 267.460, 217.890, 207.675, 230.285, 20…
    ## $ inrow_weed_biomass     <dbl> 19.000, 30.975, 0.950, 0.660, 0.495, 6.395, 1.4…
    ## $ interrow_weed_biomass  <dbl> 44.490, 0.720, 6.890, 45.735, 22.025, 19.460, 2…
    ## $ weed_biomass           <dbl> 63.490, 31.695, 7.840, 46.395, 22.520, 25.855, …
    ## $ inrow_mulch_biomass    <chr> "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",…
    ## $ interrow_mulch_biomass <chr> "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",…
    ## $ mulch_biomass          <chr> "NA", "NA", "285.95", "241.03000000000003", "30…
    ## $ bean_population        <dbl> 34.5, 39.5, 37.5, 35.0, 39.0, 33.5, 34.0, 32.0,…
    ## $ bean_yield             <dbl> 417.21, 565.54, 449.93, 412.59, 473.79, 484.04,…
    ## $ seed_weight            <dbl> 17.1200, 17.4750, 16.7525, 16.1450, 17.0475, 17…

# Initial data checks

## Site-year structure

``` r
raw_data |>
  count(location, year)
```

    ## # A tibble: 3 × 3
    ##   location  year     n
    ##   <chr>    <dbl> <int>
    ## 1 O2E       2024    20
    ## 2 O2W       2024    20
    ## 3 V         2023    20

``` r
raw_data |>
  count(location, year, treatment)
```

    ## # A tibble: 15 × 4
    ##    location  year treatment     n
    ##    <chr>    <dbl> <chr>     <int>
    ##  1 O2E       2024 RIC           4
    ##  2 O2E       2024 RIM           4
    ##  3 O2E       2024 RNO           4
    ##  4 O2E       2024 TIC           4
    ##  5 O2E       2024 TIM           4
    ##  6 O2W       2024 RIC           4
    ##  7 O2W       2024 RIM           4
    ##  8 O2W       2024 RNO           4
    ##  9 O2W       2024 TIC           4
    ## 10 O2W       2024 TIM           4
    ## 11 V         2023 RIC           4
    ## 12 V         2023 RIM           4
    ## 13 V         2023 RNO           4
    ## 14 V         2023 TIC           4
    ## 15 V         2023 TIM           4

## Identifier checks

``` r
raw_data |>
  summarise(
    n_rows = n(),
    n_unique_ids = n_distinct(id),
    n_unique_plot_records = n_distinct(location, year, block, plot),
    duplicate_ids = sum(duplicated(id)),
    duplicate_plot_records = sum(
      duplicated(paste(location, year, block, plot, sep = "_"))
    )
  )
```

    ## # A tibble: 1 × 5
    ##   n_rows n_unique_ids n_unique_plot_records duplicate_ids duplicate_plot_records
    ##    <int>        <int>                 <int>         <int>                  <int>
    ## 1     60           60                    60             0                      0

## Missing values

``` r
raw_data |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) |>
  filter(n_missing > 0)
```

    ## # A tibble: 0 × 2
    ## # ℹ 2 variables: variable <chr>, n_missing <int>

# Prepare dataset

``` r
processed_data <- raw_data |>
  mutate(
    cereal_rye_variety = recode(
      cereal_rye_variety,
      "aroostook" = "Aroostook",
      "nd gardner" = "ND Gardner"
    ),
    site_year = paste(location, year, sep = "_"),
    treatment = factor(
      treatment,
      levels = c("RNO", "RIM", "RIC", "TIM", "TIC")
    ),
    location = factor(location, levels = c("V", "O2W", "O2E")),
    year = as.integer(year),
    block = as.factor(block),
    plot = as.integer(plot)
  )
```

# Conversions

``` r
# Conversion constants
biomass_sample_area_m2 <- 0.5
row_spacing_m <- 0.762               # 30-inch rows
ha_to_ac <- 2.47105
kg_ha_to_lb_ac <- 0.892179
lb_ac_to_kg_ha <- 1.12085
lb_per_bu <- 60
g_per_lb <- 454
sqft_per_ac <- 43560
bean_yield_area_ft2 <- 16.4          # 2 m of 30-inch row

yield_moisture_factor <-
  (100 - 0.00001) / (100 - 13)

# Ensure measurement columns are numeric before calculating.
# Blank cells and common text-coded missing values become NA.
numeric_source_cols <- c(
  "bean_emergence",
  "bean_biomass",
  "inrow_weed_biomass",
  "interrow_weed_biomass",
  "weed_biomass",
  "inrow_mulch_biomass",
  "interrow_mulch_biomass",
  "mulch_biomass",
  "bean_population",
  "bean_yield",
  "seed_weight"
)

processed_data <- processed_data |>
  mutate(
    across(
      all_of(numeric_source_cols),
      ~ {
        x <- as.character(.x)
        x[x %in% c("", "NA", "N/A")] <- NA_character_
        as.numeric(x)
      }
    )
  ) |>
  mutate(
    # Biomass: all measurements converted using the 0.5 m2 basis
    bean_biomass_kg_ha =
      (bean_biomass / .env$biomass_sample_area_m2) * 10,
    bean_biomass_lb_ac =
      bean_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    weed_biomass_kg_ha =
      (weed_biomass / .env$biomass_sample_area_m2) * 10,
    weed_biomass_lb_ac =
      weed_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    inrow_weed_biomass_kg_ha =
      (inrow_weed_biomass / .env$biomass_sample_area_m2) * 10,
    inrow_weed_biomass_lb_ac =
      inrow_weed_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    interrow_weed_biomass_kg_ha =
      (interrow_weed_biomass / .env$biomass_sample_area_m2) * 10,
    interrow_weed_biomass_lb_ac =
      interrow_weed_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    mulch_biomass_kg_ha =
      (mulch_biomass / .env$biomass_sample_area_m2) * 10,
    mulch_biomass_lb_ac =
      mulch_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    inrow_mulch_biomass_kg_ha =
      (inrow_mulch_biomass / .env$biomass_sample_area_m2) * 10,
    inrow_mulch_biomass_lb_ac =
      inrow_mulch_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    interrow_mulch_biomass_kg_ha =
      (interrow_mulch_biomass / .env$biomass_sample_area_m2) * 10,
    interrow_mulch_biomass_lb_ac =
      interrow_mulch_biomass_kg_ha * .env$kg_ha_to_lb_ac,

    # Counts: plants per meter of row to plants per hectare and acre
    bean_emergence_plants_ha =
      (bean_emergence / .env$row_spacing_m) * 10000,
    bean_emergence_plants_ac =
      bean_emergence_plants_ha / .env$ha_to_ac,

    bean_population_plants_ha =
      (bean_population / .env$row_spacing_m) * 10000,
    bean_population_plants_ac =
      bean_population_plants_ha / .env$ha_to_ac,

    # Yield: adjust near-zero-moisture sample weight to 13% moisture
    bean_yield_adj_lbs_acre =
      ((bean_yield / .env$g_per_lb) /
         (.env$bean_yield_area_ft2 / .env$sqft_per_ac)) *
      .env$yield_moisture_factor,

    bean_yield_adj_bu_acre =
      bean_yield_adj_lbs_acre / .env$lb_per_bu,

    bean_yield_adj_kg_ha =
      bean_yield_adj_lbs_acre * .env$lb_ac_to_kg_ha,

    # Seed weight is already in grams per 100 seeds
    seed_weight_g_100 = seed_weight
  )
```

# Save procecessed data

``` r
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  processed_data,
  file = file.path(processed_dir, "combined_processed.rds")
)
```
