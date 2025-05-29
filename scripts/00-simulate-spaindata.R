if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(msae, tidyverse, data.table, sae, car, sf, pps)


### include a poverty line and calculate some indicators
data("incomedata")

# Draw 25% sample per province from incomedata
n_per_province <- round(0.25*(as.data.frame(table(incomedata$prov))$Freq), 0) # adjust as needed
set.seed(123)
sample_id <- stratsrs(incomedata$prov, n_per_province) 
incomedata <- incomedata[sample_id,]

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
                                 c = 0,
                                 h = 10){
  
  set.seed(seed)
  
  rho <- rho + rnorm(length(x))/h
  
  y <- c + (rho)*x 
  
  return(y)
  
  
}



incomedata <- 
  incomedata %>%
  rename(income2012 = "income") %>%
  mutate(income2013 = simulate_ar1_process(x = income2012, rho = 1.01, c = 0),
         income2014 = simulate_ar1_process(x = income2013, rho = 1.05, c = 0, h = 2))


### check correlation matrix for all numeric variables
incomedata %>%
  dplyr::select(where(is.numeric)) %>%
  cor()
### clearly not a lot of correlated variables so lets simulate some


# Create a helper function to simulate variables correlated with income
simulate_variable <- function(target, 
                              rho = 0.6, 
                              noise_sd = 1) {
  
  set.seed(123)
  # Standardize the target
  z <- scale(target)
  
  # Simulate a new variable correlated with z
  new_var <- rho * z + sqrt(1 - rho^2) * rnorm(length(z), sd = noise_sd)
  
  return(as.numeric(scale(new_var)))
  
}




### include poverty line
incomedata <-
  incomedata %>%
  mutate(povline2012 = 0.6*median(income2012),
         povline2013 = 0.6*median(income2013),
         povline2014 = 0.6*median(income2013))


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

#### create some correlated variables at the province level and simulate them back into
#### the income data

set.seed(123)

pov_dt <-
  pov_dt %>%
  mutate(
    abs = as.numeric(scale(log(poor2012 + 1) + rnorm(n(), 0, 0.07))),
    ntl = as.numeric(scale(sqrt(poor2013 + abs(min(poor2013)) + 1) + rnorm(n(), 0, 0.2))),
    aec = as.numeric(scale((as.numeric(scale(poor2014)))^2 + rnorm(n(), 0, 0.2))),
    schyrs = as.numeric(scale(qnorm((rank(poor2012) - 0.5) / n()) + rnorm(n(), 0, 0.2))),
    mkt = as.numeric(scale(poor2013 + runif(n(), -0.15, 0.15)))
  )

### merge back into the income data
incomedata <- 
  incomedata %>%
  merge(pov_dt %>%
          dplyr::select(abs, ntl, aec, schyrs, mkt, provlab, prov),
        by = c("provlab", "prov"))

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
