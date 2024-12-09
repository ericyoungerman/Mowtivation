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
cu_raw_2023 <- read_excel("~/Github/Mowtivation/raw-data/cornell_raw_2023.xlsx")
kable(head(cu_raw_2023))
```

| id | loc | year | trt | block | pllot | emerge | bbm | intrabm | interbm | totwbm | totmbm | beanden | beanyd |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | Na | 34.5 | 417.21 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | Na | 39.5 | 565.54 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 285.95 | 37.5 | 449.93 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 241.03 | 35.0 | 412.59 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 306.64 | 39.0 | 473.79 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 370.94499999999999 | 33.5 | 484.04 |

``` r
#Standardaze column names, convert to factors, check for outliers of variable**
clean_2023 <- clean_names(cu_raw_2023) |>  
  rename ('cultivation'= trt) |> 
  mutate(across(c(cultivation, block, pllot, loc), as.factor)) #|> 
  #mutate(is_outlier = totwbm < (quantile(totwbm, 0.25) - 1.5 * IQR(totwbm)) |
                       #wbm > (quantile(totwbm, 0.75) + 1.5 * IQR(totwbm)))

#select and convert data for wbm analysis
beanyd_clean_2023 <-clean_2023 |>              
  mutate(beanyd_adj_bu_acre = (((beanyd/454)/(16.4/43560))/60)* ((100-0.00001)/(100-14)))   |> 
  mutate(beanyd_adj_lbs_acre = ((beanyd/454)/(16.4/43560))* ((100-0.00001)/(100-14))) |>
  mutate(beanyd_adj_kg_ha = ((beanyd/454)/(16.4/43560))* 1.12085 *((100-0.00001)/(100-14)))
kable(head(beanyd_clean_2023)) 
```

| id | loc | year | cultivation | block | pllot | emerge | bbm | intrabm | interbm | totwbm | totmbm | beanden | beanyd | beanyd_adj_bu_acre | beanyd_adj_lbs_acre | beanyd_adj_kg_ha |
|:---|:---|---:|:---|:---|:---|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | Na | 34.5 | 417.21 | 47.30348 | 2838.209 | 3181.207 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | Na | 39.5 | 565.54 | 64.12122 | 3847.273 | 4312.216 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 285.95 | 37.5 | 449.93 | 51.01330 | 3060.798 | 3430.695 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 241.03 | 35.0 | 412.59 | 46.77967 | 2806.780 | 3145.979 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 306.64 | 39.0 | 473.79 | 53.71855 | 3223.113 | 3612.626 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 370.94499999999999 | 33.5 | 484.04 | 54.88070 | 3292.842 | 3690.782 |

<br>

## **Model testing**

### **block is fixed**

``` r
beanyd_fixed <- lm(beanyd_adj_kg_ha  ~ cultivation + block  , data = beanyd_clean_2023)

resid_panel(beanyd_fixed)
```

![](bean_yield_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### **block is random**

``` r
beanyd_ran <- lmer(beanyd_adj_kg_ha  ~ cultivation + (1|block)  , data = beanyd_clean_2023)

resid_panel(beanyd_ran)
```

![](bean_yield_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

<br>

\##**Joint test**

``` r
 beanyd_fixed |> 
  joint_tests() |> 
  kable()  
```

| model term  | df1 | df2 | F.ratio |   p.value |
|:------------|----:|----:|--------:|----------:|
| cultivation |   4 |  12 |   1.208 | 0.3575539 |
| block       |   3 |  12 |   1.291 | 0.3223078 |

<br>

# **Means comparison**

``` r
beanyd_means_2023 <- 
 emmeans(beanyd_fixed, ~  cultivation)
# Optional: Adjust for multiple comparisons (e.g., using Tukey's method)

pairwise_comparisons<- pairs(beanyd_means_2023) 
kable(head(pairwise_comparisons))
```

| contrast  |   estimate |      SE |  df |    t.ratio |   p.value |
|:----------|-----------:|--------:|----:|-----------:|----------:|
| RIC - RIM | -294.55193 | 345.403 |  12 | -0.8527775 | 0.9580254 |
| RIC - RNO |   87.40102 | 345.403 |  12 |  0.2530407 | 0.9999442 |
| RIC - TIC | -168.39709 | 345.403 |  12 | -0.4875379 | 0.9976223 |
| RIC - TIM |  406.96280 | 345.403 |  12 |  1.1782259 | 0.8378381 |
| RIM - RNO |  381.95296 | 345.403 |  12 |  1.1058181 | 0.8724211 |
| RIM - TIC |  126.15485 | 345.403 |  12 |  0.3652395 | 0.9995313 |

### **Fisher’s method for comparing means**

``` r
#mowing
cld_cultivation_fisher <-cld(emmeans(beanyd_fixed, ~  cultivation, type = "response"), Letters = letters, sort = TRUE, adjust="none", reversed=TRUE)
cld_cultivation_fisher
```

    ##  cultivation emmean  SE df lower.CL upper.CL .group
    ##  RIM           4163 244 12     3630     4695  a    
    ##  TIC           4036 244 12     3504     4569  a    
    ##  RIC           3868 244 12     3336     4400  a    
    ##  RNO           3781 244 12     3248     4313  a    
    ##  TIM           3461 244 12     2929     3993  a    
    ## 
    ## Results are averaged over the levels of: block 
    ## Confidence level used: 0.95 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

# **FIGURES**

## **Cultivation**

``` r
beanyd_clean_2023 |> 
  left_join(cld_cultivation_fisher) |> 
  ggplot(aes(x = cultivation, y = beanyd_adj_kg_ha, fill = cultivation)) +
  stat_summary(geom = "bar", fun = "mean", width = 0.7) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
  stat_summary(geom="text", fun = "MeanPlusSe", aes(label= trimws(.group)),size=6.5,vjust=-0.5) +
  labs(
    x = "Interrow weed control method",
    y = expression("Soybean yield" ~ (kg ~ ha^{-1})),
    title = str_c("The influence of interrow weed control method on soybean yield"),
    subtitle = expression(italic("P = 0.01"))) +
  
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
ggsave("beanyd_plot_cutivation.png", width = 8, height = 6, dpi = 300)
```
