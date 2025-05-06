################################################################################
############### PREPARE SPAIN DATA FOR MFH SPATIO-TEMPORAL SAE #################
################################################################################

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(msae, tidyverse, data.table, sae, car, sf)


### include a poverty line and calculate some indicators
data("incomedata")



simulate_correlated_vector <- function(x, 
                                       rho){
  
  set.seed(123)
  
  xstd <- scale(x)[,1]
  
  ## retransform since we are going to be taking the absolute value
  ## of the income variable
  rho <- rho^2 
  
  # y <- rho * xstd + sqrt(1 - rho^2) * abs(rnorm(length(x)))
  # 
  # y <- y * sd(x) + mean(y)
  
  
  
  return(y)
    
}

simulate_ar1_process <- function(x, 
                                 rho,
                                 seed = 123,
                                 c = 0){
  
  set.seed(seed)
  
  rho <- rho * rnorm(length(x))/10
  
  y <- c + (rho + 1)*x 
  
  return(y)
  
  
}



incomedata <- 
  incomedata %>%
  rename(income2012 = "income") %>%
  mutate(income2013 = simulate_ar1_process(x = income2012, rho = 0.8, c = 50),
         income2014 = simulate_ar1_process(x = income2013, rho = 0.7, c = 0))


#### include some variables

set.seed(123)

incomedata <- 
  incomedata %>%
  mutate(
    ntl = 0.7 * income2012 + 0.3 * income2013*rnorm(n = nrow(.)), 
    yos = 0.6 * income2013 + 0.4 * income2014*rnorm(n = nrow(.)), 
    abc = 0.8 * income2012 + 0.2 * income2014*rnorm(n = nrow(.)),
    poi = 0.8 * income2012 + 0.7 * income2012*rnorm(n = nrow(.)),
    ips = income2014 + 0.9 * income2013*rnorm(n = nrow(.)),
    sam = 0.4 * income2013 + 0.6 * income2012*rnorm(n = nrow(.))
  )
### include poverty line
incomedata <- 
  incomedata %>%
  mutate(povline2012 = 6557.143,
         povline2013 = 6557.143 * 1.1,
         povline2014 = 6557.143 * 1.2)


### lets test the correlation in poverty rates across the years
pov_dt <- 
incomedata %>%
  mutate(poor2012 = ifelse(income2012 < povline2012, 1, 0),
         poor2013 = ifelse(income2013 < povline2013, 1, 0),
         poor2014 = ifelse(income2014 < povline2014, 1, 0)) %>%
  as.data.table() %>%
  .[, lapply(.SD, weighted.mean, w = weight), 
    .SDcols = c("poor2012", "poor2013", "poor2014"),
    by = c("provlab", "prov")]
  


### combine the datasets and include a shapefile for proximity estimation
shp_dt <- sf::read_sf("data/shapes/georef-spain-provincia-millesime.shp")

spain_dt <- 
  incomedata[, c("prov", "provlab")] %>% 
  unique() %>%
  merge(shp_dt %>% 
          mutate(prov_code = as.integer(prov_code)) %>% 
          dplyr::select(prov_code, prov_name) %>% 
          unique(.),
        by.x = "prov",
        by.y = "prov_code")

saveRDS(spain_dt, "data/shapes/spainshape.RDS")



saveRDS(incomedata, "data/incomedata.RDS")















