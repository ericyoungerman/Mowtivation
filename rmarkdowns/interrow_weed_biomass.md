interrow weed biomass
================

# Load libraries

``` r
#Set work directory
setwd("/Users/ey239/Github/Mowtivation/rmarkdowns")

#Load packages 
library(tidyverse) ##install.packages("tidyverse")
library(knitr)
library(patchwork) ##install.packages("patchwork")
library(skimr)     ##install.packages("skimr")
library(readxl)
library(janitor) ##install.packages("janitor")

library(kableExtra) ##install.packages("kableExtra")
library(webshot) ##install.packages("webshot")
webshot::install_phantomjs()
library(viridis) ##install.packages("viridis")
library(lme4) ##install.packages("lme4")
library(lmerTest) ##install.packages("lmerTest")
library(emmeans) ##install.packages("emmeans")
library(rstatix) ##install.packages("rstatix")
#library(Matrix) ##install.packages("Matrix")
library(multcomp) ##install.packages("multcomp")
library(multcompView) ##install.packages("multcompView")
library(ggResidpanel) ##install.packages("ggResidpanel")
#library(car)
#library(TMB)  ##install.packages("TMB")
#library(glmmTMB)  ##install.packages("glmmTMB")
#library(DHARMa)  ##install.packages("DHARMa")

#Load Functions
MeanPlusSe<-function(x) mean(x)+plotrix::std.error(x)

find_logw0=function(x){c=trunc(log(min(x[x>0],na.rm=T)))
d=exp(c)
return(d)}
```

<br>

# Load and Clean Data

## Load data

``` r
combined_raw <- read_excel("~/Github/Mowtivation/raw-data/All Treatments/combined_raw.xlsx")
kable(head(combined_raw))
```

| id | location | year | treatment | block | plot | bean_emergence | bean_biomass | intrarow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 37.5 | 449.93 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35 | 412.59 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39 | 473.79 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 |

## Clean data

``` r
#Standardaze column names, convert to factors, check for outliers of variable**
clean_combined <- clean_names(combined_raw) |>  
  rename ('weed_control'= treatment) |> 
  mutate(across(c(weed_control, block, plot, location, year), as.factor)) #|> 
  #mutate(is_outlier = totwbm < (quantile(totwbm, 0.25) - 1.5 * IQR(totwbm)) |
                       #wbm > (quantile(totwbm, 0.75) + 1.5 * IQR(totwbm)))

#select and convert data for wbm analysis
interrow_weed_biomass_clean <-clean_combined |>              
  mutate(intrarow_weed_biomass_grams_meter = (interrow_weed_biomass * 2)) |> 
  mutate(interrow_weed_biomass_kg_ha = ((interrow_weed_biomass/0.5) *(10000))/(1000)) |>
  mutate(interrow_weed_biomass_lbs_ac = (((interrow_weed_biomass/0.5) *(10000))/(1000))* 0.892179)
kable(head(interrow_weed_biomass_clean)) 
```

| id | location | year | weed_control | block | plot | bean_emergence | bean_biomass | intrarow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | intrarow_weed_biomass_grams_meter | interrow_weed_biomass_kg_ha | interrow_weed_biomass_lbs_ac |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|:---|:---|---:|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 88.98 | 889.8 | 793.86087 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 1.44 | 14.4 | 12.84738 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 37.5 | 449.93 | 13.78 | 137.8 | 122.94227 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35 | 412.59 | 91.47 | 914.7 | 816.07613 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39 | 473.79 | 44.05 | 440.5 | 393.00485 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 38.92 | 389.2 | 347.23607 |

# Model testing

Block is random Tyler is under the impression that block should always
be random and that post-hoc comparisons should use TUKEY rather the
Fischer. Fisher is bogus apparently.

## Lmer

``` r
random <- lmer(interrow_weed_biomass_lbs_ac  ~ location+weed_control + location:weed_control +(1|location:block) , data = interrow_weed_biomass_clean)

resid_panel(random)
```

![](interrow_weed_biomass_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

\##Joint test (anova)

``` r
random |> 
  joint_tests() |> 
  kable()  
```

|     | model term            | df1 | df2 | F.ratio |   p.value |
|:----|:----------------------|----:|----:|--------:|----------:|
| 1   | location              |   2 |   9 |   3.924 | 0.0595212 |
| 3   | weed_control          |   4 |  36 |   4.136 | 0.0073744 |
| 2   | location:weed_control |   8 |  36 |   1.079 | 0.3996164 |

<br>

## Means comparison

### Weed-control (S)

``` r
means_weed_control <- emmeans(random, ~  weed_control)
pairwise_comparisons_weed_control<- pairs(means_weed_control) 
kable(head(pairwise_comparisons_weed_control))
```

| contrast  |   estimate |       SE |  df |    t.ratio |   p.value |
|:----------|-----------:|---------:|----:|-----------:|----------:|
| RIC - RIM |   98.38504 | 69.04058 |  36 |  1.4250321 | 0.6555723 |
| RIC - RNO |  -68.77213 | 69.04058 |  36 | -0.9961117 | 0.9061227 |
| RIC - TIC |  141.23194 | 69.04058 |  36 |  2.0456366 | 0.2562870 |
| RIC - TIM |  -81.66412 | 69.04058 |  36 | -1.1828423 | 0.8142297 |
| RIM - RNO | -167.15717 | 69.04058 |  36 | -2.4211438 | 0.1176228 |
| RIM - TIC |   42.84690 | 69.04058 |  36 |  0.6206045 | 0.9903726 |

br\>

### Location (NS)

``` r
means_location <- emmeans(random, ~  location)
pairwise_comparisons_location<- pairs(means_location) 
kable(head(pairwise_comparisons_location))
```

| contrast | estimate | SE | df | t.ratio | p.value |
|:---|---:|---:|---:|---:|---:|
| field O2 east - field O2 west | -54.88239 | 68.73638 | 9 | -0.7984475 | 0.8292015 |
| field O2 east - field x | -187.28175 | 68.73638 | 9 | -2.7246379 | 0.0686539 |
| field O2 west - field x | -132.39936 | 68.73638 | 9 | -1.9261904 | 0.2369423 |
| \### Location | weed-control (S) |  |  |  |  |

``` r
means_weed_control_location <- emmeans(random, ~  weed_control|location)
pairwise_comparisons_weed_control_location<- pairs(means_weed_control_location) 
kable(head(pairwise_comparisons_weed_control_location))
```

| contrast  | location      |    estimate |       SE |  df |    t.ratio |   p.value |
|:----------|:--------------|------------:|---------:|----:|-----------:|----------:|
| RIC - RIM | field O2 east |  -19.962505 | 119.5818 |  36 | -0.1669360 | 0.9999948 |
| RIC - RNO | field O2 east | -109.760322 | 119.5818 |  36 | -0.9178682 | 0.9343134 |
| RIC - TIC | field O2 east |  -25.560928 | 119.5818 |  36 | -0.2137527 | 0.9999775 |
| RIC - TIM | field O2 east | -101.953755 | 119.5818 |  36 | -0.8525859 | 0.9531223 |
| RIM - RNO | field O2 east |  -89.797816 | 119.5818 |  36 | -0.7509322 | 0.9745292 |
| RIM - TIC | field O2 east |   -5.598423 | 119.5818 |  36 | -0.0468167 | 1.0000000 |

## Tukey compact letter display

### Weed-control (S)

``` r
#weed control
cld_weed_control_tukey <-cld(emmeans(random, ~  weed_control , type = "response"), Letters = letters, sort = TRUE, reversed=TRUE)
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
cld_weed_control_tukey
```

    ##  weed_control emmean   SE   df lower.CL upper.CL .group
    ##  TIM           233.2 51.9 42.7    128.5      338  a    
    ##  RNO           220.3 51.9 42.7    115.6      325  a    
    ##  RIC           151.5 51.9 42.7     46.8      256  ab   
    ##  RIM            53.2 51.9 42.7    -51.5      158  ab   
    ##  TIC            10.3 51.9 42.7    -94.4      115   b   
    ## 
    ## Results are averaged over the levels of: location 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

# Figures

## Weed control

``` r
interrow_weed_biomass_clean |> 
  left_join(cld_weed_control_tukey) |> 
  ggplot(aes(x = weed_control, y = interrow_weed_biomass_lbs_ac, fill = weed_control)) +
  stat_summary(geom = "bar", fun = "mean", width = 0.7) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
  stat_summary(geom="text", fun = "MeanPlusSe", aes(label= trimws(.group)),size=6.5,vjust=-0.5) +
  labs(
    x = "Interrow weed control",
    y = expression("Weed biomass" ~ (lb ~ A^{-1})),
    title = str_c("Influence of interrow weed control on interrow weed biomass"),
    subtitle = expression(italic("P = 0.007"))) +
  
  scale_x_discrete(labels = c("Rolled,\nhigh-residue\ncultivation",
                              "Rolled,\ninterrow\nmowing",
                              "Rolled,\nno additional\nweed control",
                          "Tilled,\nstandard\ncultivation",
                              "Tilled,\ninterrow\nmowing")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_viridis(discrete = TRUE, option = "D", direction = -1, end = 0.9, begin = 0.1) +
   theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )
```

![](interrow_weed_biomass_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave("interrow_weed_biomass_weed_control_lbA.png", width = 8, height = 6, dpi = 300)
```
