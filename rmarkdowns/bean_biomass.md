bean biomass
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
  #mutate(is_outlier = bbm < (quantile(bbm, 0.25) - 1.5 * IQR(bbm)) |
                       #bbm > (quantile(bbm, 0.75) + 1.5 * IQR(bbm)))

#select and convert data for wbm analysis
  
bbm_clean_2023 <-clean_2023 |>              
  
  mutate(bbm_grams_meter = (bbm * 2)) |> 
  mutate(bbm_kg_ha = ((bbm/0.5) *(10000))/(1000)) |>
  mutate(bbm_lbs_ac = (((bbm/0.5) *(10000))/(1000))* 0.892179)
kable(head(bbm_clean_2023)) 
```

| id | loc | year | cultivation | block | pllot | emerge | bbm | intrabm | interbm | totwbm | totmbm | beanden | beanyd | bbm_grams_meter | bbm_kg_ha | bbm_lbs_ac |
|:---|:---|---:|:---|:---|:---|---:|---:|---:|---:|---:|:---|---:|---:|---:|---:|---:|
| CU_B1_P101 | field x | 2023 | TIM | 1 | 101 | 46.5 | 223.740 | 19.000 | 44.490 | 63.490 | Na | 34.5 | 417.21 | 447.48 | 4474.8 | 3992.323 |
| CU_B1_P102 | field x | 2023 | TIC | 1 | 102 | 42.5 | 267.460 | 30.975 | 0.720 | 31.695 | Na | 39.5 | 565.54 | 534.92 | 5349.2 | 4772.444 |
| CU_B1_P103 | field x | 2023 | RIM | 1 | 103 | 36.5 | 217.890 | 0.950 | 6.890 | 3.920 | 285.95 | 37.5 | 449.93 | 435.78 | 4357.8 | 3887.938 |
| CU_B1_P104 | field x | 2023 | RNO | 1 | 104 | 41.0 | 207.675 | 0.660 | 45.735 | 46.395 | 241.03 | 35.0 | 412.59 | 415.35 | 4153.5 | 3705.665 |
| CU_B1_P105 | field x | 2023 | RIC | 1 | 105 | 41.0 | 230.285 | 0.495 | 22.025 | 22.520 | 306.64 | 39.0 | 473.79 | 460.57 | 4605.7 | 4109.109 |
| CU_B1_P201 | field x | 2023 | RIC | 2 | 201 | 36.5 | 208.105 | 6.395 | 19.460 | 25.855 | 370.94499999999999 | 33.5 | 484.04 | 416.21 | 4162.1 | 3713.338 |
