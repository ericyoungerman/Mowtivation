Soybean emergence
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

<br> \# Load and clean data

\##Load data

``` r
combined_raw <- read_excel("~/Github/Mowtivation/raw-data/All Treatments/combined_raw.xlsx")
kable(head(combined_raw))
```

| id | location | year | treatment | block | plot | bean_emergence | bean_biomass | intrarow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | seed_weight |
|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 17.119999999999997 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 17.475000000000001 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 37.5 | 449.93 | 16.752499999999998 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35 | 412.59 | 16.145 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39 | 473.79 | 17.047499999999999 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 17.149999999999999 |

## Clean data

``` r
#Standardaze column names, convert to factors, check for outliers of variable**
clean_combined <- clean_names(combined_raw) |>  
  rename ('weed_control'= treatment) |> 
  mutate(across(c(weed_control, block, plot, location, year), as.factor)) #|> 
  #mutate(is_outlier = totwbm < (quantile(totwbm, 0.25) - 1.5 * IQR(totwbm)) |
                       #wbm > (quantile(totwbm, 0.75) + 1.5 * IQR(totwbm)))

#select and convert data for wbm analysis
bean_emergence_clean <-clean_combined |>              
  mutate(bean_emergence_two_meter = (bean_emergence * 2)) |> 
  mutate(bean_emergence_acre = (((bean_emergence/0.762) *10000)/2.471))

kable(head(bean_emergence_clean)) 
```

| id | location | year | weed_control | block | plot | bean_emergence | bean_biomass | intrarow_weed_biomass | interrow_weed_biomass | weed_biomass | bean_population | bean_yield | seed_weight | bean_emergence_two_meter | bean_emergence_acre |
|:---|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|:---|:---|:---|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | 34.5 | 417.21 | 17.119999999999997 | 93 | 246959.2 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | 39.5 | 565.54 | 17.475000000000001 | 85 | 225715.4 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 37.5 | 449.93 | 16.752499999999998 | 73 | 193849.7 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 35 | 412.59 | 16.145 | 82 | 217749.0 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 39 | 473.79 | 17.047499999999999 | 82 | 217749.0 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 33.5 | 484.04 | 17.149999999999999 | 73 | 193849.7 |

\#Model testing

``` r
random <- lmer(bean_emergence_acre  ~ location+weed_control + location:weed_control +(1|location:block) , data = bean_emergence_clean)

resid_panel(random)
```

![](bean_emergence_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

\##Joint test (anova)

``` r
random |> 
  joint_tests() |> 
  kable()  
```

|     | model term            | df1 | df2 | F.ratio |   p.value |
|:----|:----------------------|----:|----:|--------:|----------:|
| 1   | location              |   2 |   9 |   7.721 | 0.0111542 |
| 3   | weed_control          |   4 |  36 |   4.487 | 0.0048141 |
| 2   | location:weed_control |   8 |  36 |   1.218 | 0.3165498 |

## Means comparison

### Weed-control (S)

However, since the interrow weed treatments were implemented after
emergence. Significance is due to pre-emergent weed control brought
about by rolling and no rolling

``` r
means_weed_control <- 
 emmeans(random, ~  weed_control)

pairwise_comparisons_weed_control<- pairs(means_weed_control) 
kable(head(pairwise_comparisons_weed_control))
```

| contrast  |   estimate |       SE |  df |    t.ratio |   p.value |
|:----------|-----------:|---------:|----:|-----------:|----------:|
| RIC - RIM |  -2434.186 | 7484.258 |  36 | -0.3252408 | 0.9997370 |
| RIC - RNO |   6638.689 | 7484.258 |  36 |  0.8870203 | 0.9437231 |
| RIC - TIC | -18145.749 | 7484.258 |  36 | -2.4245221 | 0.1167376 |
| RIC - TIM | -18145.749 | 7484.258 |  36 | -2.4245221 | 0.1167376 |
| RIM - RNO |   9072.874 | 7484.258 |  36 |  1.2122611 | 0.7968936 |
| RIM - TIC | -15711.563 | 7484.258 |  36 | -2.0992813 | 0.2311649 |

### Location (S)

``` r
means_location <- 
 emmeans(random, ~  location)

pairwise_comparisons_location<- pairs(means_location) 
kable(head(pairwise_comparisons_location))
```

| contrast                      |  estimate |       SE |  df |    t.ratio |   p.value |
|:------------------------------|----------:|---------:|----:|-----------:|----------:|
| field O2 east - field O2 west | -4115.987 | 6502.753 |   9 | -0.6329607 | 0.9042459 |
| field O2 east - field x       | 19783.292 | 6502.753 |   9 |  3.0422950 | 0.0413154 |
| field O2 west - field x       | 23899.279 | 6502.753 |   9 |  3.6752557 | 0.0152606 |

## Tukey compact letter display

### Weed-control (S)

Significant but related to rolled vs. tilled rather than interrow weed
control

``` r
cld_weed_control_tukey <-cld(emmeans(random, ~  weed_control, type = "response"), Letters = letters, sort = TRUE, reversed=TRUE)
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
cld_weed_control_tukey
```

    ##  weed_control emmean   SE   df lower.CL upper.CL .group
    ##  TIM          243640 5430 44.6   232706   254574  a    
    ##  TIC          243640 5430 44.6   232706   254574  a    
    ##  RIM          227928 5430 44.6   216995   238862  ab   
    ##  RIC          225494 5430 44.6   214560   236428  ab   
    ##  RNO          218855 5430 44.6   207922   229789   b   
    ## 
    ## Results are averaged over the levels of: location 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

### Location (S)

``` r
#location
cld_location_tukey <-cld(emmeans(random, ~  location , type = "response"), Letters = letters, sort = TRUE, reversed=TRUE)
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
cld_location_tukey
```

    ##  location      emmean   SE df lower.CL upper.CL .group
    ##  field O2 west 241250 4600  9   230848   251652  a    
    ##  field O2 east 237134 4600  9   226732   247536  a    
    ##  field x       217351 4600  9   206949   227752   b   
    ## 
    ## Results are averaged over the levels of: weed_control 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

### Weed-control:Location (NS)

``` r
cld_weed_control_location_tukey <-cld(emmeans(random, ~  weed_control|location, type = "response"), Letters = letters, sort = TRUE, reversed=TRUE)
cld_weed_control_location_tukey
```

    ## location = field O2 east:
    ##  weed_control emmean   SE   df lower.CL upper.CL .group
    ##  TIC          258909 9400 44.6   239971   277846  a    
    ##  TIM          238329 9400 44.6   219391   257267  a    
    ##  RIM          237001 9400 44.6   218064   255939  a    
    ##  RNO          227043 9400 44.6   208106   245981  a    
    ##  RIC          224388 9400 44.6   205450   243325  a    
    ## 
    ## location = field O2 west:
    ##  weed_control emmean   SE   df lower.CL upper.CL .group
    ##  TIM          253598 9400 44.6   234660   272536  a    
    ##  RIM          244968 9400 44.6   226030   263905  a    
    ##  RIC          242976 9400 44.6   224038   261914  a    
    ##  TIC          240321 9400 44.6   221383   259258  a    
    ##  RNO          224388 9400 44.6   205450   243325  a    
    ## 
    ## location = field x:
    ##  weed_control emmean   SE   df lower.CL upper.CL .group
    ##  TIM          238993 9400 44.6   220055   257930  a    
    ##  TIC          231690 9400 44.6   212753   250628  a    
    ##  RIC          209119 9400 44.6   190181   228056  a    
    ##  RNO          205135 9400 44.6   186198   224073  a    
    ##  RIM          201816 9400 44.6   182879   220754  a    
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same.

\#Figures \## Location (S)

``` r
bean_emergence_clean |> 
  left_join(cld_location_tukey) |> 
  ggplot(aes(x = location, y = bean_emergence_acre, fill = location)) +
  stat_summary(geom = "bar", fun = "mean", width = 0.7) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
  stat_summary(geom="text", fun = "MeanPlusSe", aes(label= trimws(.group)),size=6.5,vjust=-0.5) +
  
  labs(
    x = "Location",
    y = expression("Soybean emergence" ~ (plants * "/" * a)),
    title = str_c("The influence of location on soybean emergence"),
    subtitle = expression(italic("P = 0.01"))) +
   scale_x_discrete(labels = c("Field O2 East ",
                              "Field O2 West",
                              "Field X")) +
  scale_y_continuous(labels = scales::label_comma(),expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_viridis(discrete = TRUE, option = "D", direction = -1, end = 0.9, begin = 0.1) +
   theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )
```

![](bean_emergence_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave("bean_emergence_location_plantsAc.png", width = 8, height = 6, dpi = 300)
```
