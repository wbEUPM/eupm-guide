################################################################################
############### PREPARE SPAIN DATA FOR MFH SPATIO-TEMPORAL SAE #################
################################################################################

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(msae, tidyverse, data.table, sae, car, sf)


### include a poverty line and calculate some indicators
data("incomedata")


# Function to simulate data for a new year
simulate_year <- function(df, 
                          income_growth_mean = 0.05, 
                          income_growth_sd = 0.02, 
                          age_increment = 1) {
  new_df <- df

  # Simulate income increase: apply household-specific growth
  growth_factors <- rnorm(nrow(df), mean = income_growth_mean, sd = income_growth_sd)
  new_df$income <- round(df$income * (1 + growth_factors), 2)
  
  # Gender remains the same
  return(new_df)
  
}



incomedata <- 
  incomedata %>%
  rename(income2012 = "income") %>%
  cbind(simulate_year(df = incomedata) %>%
          rename(income2013 = "income") %>%
          dplyr::select(income2013)) %>%
  cbind(simulate_year(incomedata, 
                      income_growth_mean = 0.1, 
                      income_growth_sd = 0.05) %>%
          rename(income2014 = "income") %>%
          dplyr::select(income2014))


### include poverty line
incomedata <- 
  incomedata %>%
  mutate(povline2012 = 6557.143,
         povline2013 = 6557.143 * 1.01,
         povline2014 = 6557.143 * 1.02)


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















