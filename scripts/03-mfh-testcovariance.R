################################################################################
################ TEST COVARIANCE MATRIX TYPES FOR THE MFH MODELS ###############
################################################################################

## read in the data
if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(sf, data.table, tidyverse, car, msae, 
               sae, survey, spdep, knitr, MASS, caret,
               purrr, here)

# Loading locally-developed
list.files("R", pattern = "*.R$", full.names = TRUE, ignore.case = TRUE) |>
  walk(~ suppressMessages(source(.x)))


income_dt <- readRDS(here::here("data/incomedata_survey.RDS"))

prov_dt <- readRDS(here::here("data/shapes/simadmin.RDS"))

shp_dt <- readRDS(here::here("data/shapes/spainshape.RDS"))

avg_psu_size <- 12

income_dt <- 
  income_dt %>%
  group_by(provlab, prov) %>%
  group_modify(~ {
    n <- nrow(.x)
    n_psus <- max(2, ceiling(n / avg_psu_size))
    .x$psu <- sample(1:n_psus, size = n, replace = TRUE)
    .x
  }) %>%
  ungroup()






















