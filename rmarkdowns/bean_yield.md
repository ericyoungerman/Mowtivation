Bean Yield
================

# **Load libraries**

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

# **Load and Clean Data**

### **Load individual datasets**

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

``` r
#Standardaze column names, convert to factors, check for outliers of variable**
clean_combined <- clean_names(combined_raw) |>  
  rename ('weed_control'= treatment) |> 
  mutate(across(c(weed_control, block, plot, location, year), as.factor)) #|> 
  #mutate(is_outlier = totwbm < (quantile(totwbm, 0.25) - 1.5 * IQR(totwbm)) |
                       #wbm > (quantile(totwbm, 0.75) + 1.5 * IQR(totwbm)))

#select and convert data for wbm analysis
  bean_yield_clean <- clean_combined |>  
  mutate(bean_yield = as.numeric(bean_yield)) |>  # Convert beanyd to numeric
  filter(!is.na(bean_yield)) |>  # Exclude rows with NA in beanyd
  mutate(
    bean_yield_adj_bu_acre = (((bean_yield / 454) / (16.4 / 43560)) / 60) * ((100 - 0.00001) / (100 - 14)),
    bean_yield_adj_lbs_acre = ((bean_yield / 454) / (16.4 / 43560)) * ((100 - 0.00001) / (100 - 14)),
    bean_yield_adj_kg_ha = ((bean_yield / 454) / (16.4 / 43560)) * 1.12085 * ((100 - 0.00001) / (100 - 14))
  )
```

    ## Warning: There was 1 warning in `mutate()`.
    ## ℹ In argument: `bean_yield = as.numeric(bean_yield)`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

``` r
kable(head(bean_yield_clean)) 
```

| id | location | year | weed_control | block | plot | bean_emergence | bean_biomass | intrarow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | bean_yield_adj_bu_acre | bean_yield_adj_lbs_acre | bean_yield_adj_kg_ha |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 47.30348 | 2838.209 | 3181.207 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 64.12122 | 3847.273 | 4312.216 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 37.5 | 449.93 | 51.01330 | 3060.798 | 3430.695 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35 | 412.59 | 46.77967 | 2806.780 | 3145.979 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39 | 473.79 | 53.71855 | 3223.113 | 3612.626 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 54.88070 | 3292.842 | 3690.782 |

<br>

## **Model testing**

### **block is fixed**

\#Ask tyler about model format, should block always be fixed, etc.
should location be nested in year?

``` r
fixed <- lmer(bean_yield_adj_kg_ha ~ year*weed_control*location + block 
                     +(1 | year/location), data = bean_yield_clean)
```

    ## fixed-effect model matrix is rank deficient so dropping 15 columns / coefficients

    ## Warning in as_lmerModLT(model, devfun): Model may not have converged with 2
    ## eigenvalues close to zero: 2.2e-09 8.1e-10

``` r
resid_panel(fixed)
```

![](bean_yield_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### **block is random**

``` r
random <- lmer(bean_yield_adj_kg_ha  ~ weed_control + (1|block)  , data = bean_yield_clean)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
resid_panel(random)
```

![](bean_yield_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

<br>

\##**Joint test**

``` r
 fixed |> 
  joint_tests() |> 
  kable()  
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     location %in% year

|     | model term                 | df1 |      df2 | F.ratio |   p.value | note |
|:----|:---------------------------|----:|---------:|--------:|----------:|:-----|
| 1   | year                       |   1 | 90870.37 |   1.452 | 0.2282015 |      |
| 5   | weed_control               |   4 |    41.00 |   1.335 | 0.2734328 |      |
| 6   | block                      |   3 |    41.00 |   0.123 | 0.9458155 |      |
| 2   | year:weed_control          |   4 |    41.00 |   0.652 | 0.6290763 |      |
| 4   | year:location              |   1 | 16912.91 |   0.074 | 0.7852912 | e    |
| 3   | year:weed_control:location |   4 |    41.00 |   1.679 | 0.1733297 | e    |

<br>

# **Means comparison**

``` r
means <- emmeans(fixed, ~  weed_control)
# Optional: Adjust for multiple comparisons (e.g., using Tukey's method)

pairwise_comparisons<- pairs(means) 
kable(head(pairwise_comparisons))
```

| contrast  |   estimate |       SE |  df |    t.ratio |   p.value |
|:----------|-----------:|---------:|----:|-----------:|----------:|
| RIC - RIM |  -81.02477 | 206.4819 |  41 | -0.3924062 | 0.9992229 |
| RIC - RNO |   96.12206 | 203.4762 |  41 |  0.4723996 | 0.9977921 |
| RIC - TIC | -127.08890 | 203.4762 |  41 | -0.6245887 | 0.9899819 |
| RIC - TIM |  293.07460 | 203.4762 |  41 |  1.4403388 | 0.6420372 |
| RIM - RNO |  177.14683 | 206.4819 |  41 |  0.8579292 | 0.9514074 |
| RIM - TIC |  -46.06413 | 206.4819 |  41 | -0.2230904 | 0.9999709 |

### **Fisher’s method for comparing means**

``` r
#weed_control
cld_weed_control_fisher <-cld(emmeans(fixed, ~  weed_control, type = "response"), Letters = letters, sort = TRUE, adjust="none", reversed=TRUE)
```

    ## NOTE: A nesting structure was detected in the fitted model:
    ##     location %in% year

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
cld_weed_control_fisher
```

    ##  weed_control emmean  SE   df lower.CL upper.CL .group
    ##  TIC            4560 462 4371     3654     5467  a    
    ##  RIM            4514 464 3939     3605     5423  ab   
    ##  RIC            4433 462 4371     3527     5339  ab   
    ##  RNO            4337 462 4371     3431     5243  ab   
    ##  TIM            4140 462 4371     3234     5046   b   
    ## 
    ## Results are averaged over the levels of: block, location, year 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

# **FIGURES**

## **Cultivation**

``` r
bean_yield_clean |> 
  left_join(cld_weed_control_fisher) |> 
  ggplot(aes(x = weed_control, y = bean_yield_adj_kg_ha, fill = weed_control)) +
  stat_summary(geom = "bar", fun = "mean", width = 0.7) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
  #stat_summary(geom="text", fun = "MeanPlusSe", aes(label= trimws(.group)),size=6.5,vjust=-0.5) +
  labs(
    x = "Method of interrow weed control",
    y = expression("Soybean yield" ~ (kg ~ ha^{-1})),
    title = str_c("The influence of the method of interrow weed control on soybean yield"),
    subtitle = expression(italic("Not signficant"))) +
  
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

![](bean_yield_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave("beanydall_plot_cutivation.png", width = 8, height = 6, dpi = 300)
```
