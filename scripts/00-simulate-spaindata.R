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

### compute the variable set
set.seed(123)

dt <- 
incomedata |>
  group_by(provlab) |>
  sample_frac(size = 0.25)

survey_dt <- 
  dt |>
  mutate(weight = weight * 10) |>
  ungroup() |>
  dplyr::select(provlab, prov, starts_with("income"), starts_with("povline"), weight)


#' Create Interaction Terms Between a Variable and a List of Other Variables
#'
#' This function generates interaction terms between a specified variable (`interacter_var`) 
#' and a list of other variables (`var_list`) in a given dataset. The resulting interaction 
#' terms are returned in a new data frame.
#'
#' @param dt A `data.frame` or `data.table` containing the data.
#' @param interacter_var A string specifying the name of the variable to interact with each variable in `var_list`.
#' @param var_list A character vector of variable names in `dt` to be interacted with `interacter_var`.
#'
#' @return A `data.frame` containing the interaction terms. Each column is named using the pattern `var_X_interacter_var`.
#'
#' @examples
#' \dontrun{
#' dt <- data.frame(a = 1:5, b = 6:10, c = 11:15)
#' create_interactions(dt, interacter_var = "a", var_list = c("b", "c"))
#' }
#'
#' @export


create_interactions <- function(dt, interacter_var, var_list) {
  # Ensure dt is a data.table
  if (!"data.frame" %in% class(dt)) {
    dt <- as.data.table(dt)
  }
  
  # Check if interacter_var exists in the dataset
  if (!(interacter_var %in% names(dt))) {
    stop(paste(interacter_var, "not found in dataset"))
  }
  
  # Check if var_list contains valid variables that exist in the dataset
  if (any(!var_list %in% names(dt))) {
    stop("Some variables in var_list are not found in the dataset.")
  }
  
  # Create an empty data.table to store interactions
  int_dt <- data.frame(matrix(nrow = nrow(dt)))
  
  # Loop over var_list to create interaction terms
  for (var in var_list) {
    interaction_name <- paste0(var, "_X_", interacter_var)
    int_dt[[interaction_name]] <- dt[[var]] * dt[[interacter_var]]
  }
  
  int_dt <- int_dt[, -1]
  
  return(int_dt)
}

### applying the create_interactions() functions
incomedata <- 
  incomedata %>%
  cbind(create_interactions(dt = .,
                            interacter_var = "gen",
                            var_list = c("age2", "age3", "age4", "age5",
                                         "educ1", "educ2", "educ3", "nat1",
                                         "labor1", "labor2", "labor3"))) %>%
  cbind(create_interactions(dt = .,
                            interacter_var = "educ3",
                            var_list = c("age2", "age3", "age4", "age5",
                                         "nat1", "labor1", "labor2", "labor3"))) 

candidate_vars <- colnames(incomedata)[!colnames(incomedata) %in% 
                                         c("provlab", "prov", "income",
                                           "weight", "povline", "y",
                                           "poverty", "y0", "y1", "y2",
                                           "p0_prov", "p1_prov", "p2_prov",
                                           "ac", "nat", "educ", "labor",
                                           "age")]

candidate_vars <- candidate_vars[!grepl("sampsize|prov|poverty|income|povline|^v[0-9]|^y[0-9]", 
                                        candidate_vars)]

### compute the variable set
set.seed(123)

prov_dt <- 
incomedata |>
  group_by(provlab) |>
  sample_frac(size = 0.25) |>
  mutate(weight = weight * 10) |>
  ungroup() |>
  group_by(prov, provlab) |>
  summarize(
    across(
      any_of(candidate_vars),
      ~ weighted.mean(x = ., w = weight, na.rm = TRUE),
      .names = "{.col}"
    )
  ) 


income_dt <- readRDS("data/incomedata_sample.RDS")

income_dt <- 
  income_dt %>%
  cbind(create_interactions(dt = .,
                            interacter_var = "gen",
                            var_list = c("age2", "age3", "age4", "age5",
                                         "educ1", "educ2", "educ3", "nat1",
                                         "labor1", "labor2", "labor3"))) %>%
  cbind(create_interactions(dt = .,
                            interacter_var = "educ3",
                            var_list = c("age2", "age3", "age4", "age5",
                                         "nat1", "labor1", "labor2", "labor3"))) 



## create the candidate variables
candidate_vars <- colnames(income_dt)[!colnames(income_dt) %in% 
                                        c("provlab", "prov", "income",
                                          "weight", "povline", "y",
                                          "poverty", "y0", "y1", "y2",
                                          "p0_prov", "p1_prov", "p2_prov",
                                          "ac", "nat", "educ", "labor",
                                          "age")]

candidate_vars <- candidate_vars[!grepl("sampsize|prov|poverty|income|povline|^v[0-9]|^y[0-9]", 
                                        candidate_vars)]

### computing the province level data
prov_dt <- 
  income_dt |>
  group_by(prov, provlab) |>
  summarize(
    across(
      any_of(candidate_vars),
      ~ weighted.mean(x = ., w = weight, na.rm = TRUE),
      .names = "{.col}"
    )
  ) 

### select the poverty data
survey_dt <- 
 income_dt |>
  dplyr::select(provlab, prov, starts_with("income"), starts_with("povline"), weight)

### create the set of model variables





saveRDS(survey_dt, "data/incomedata_survey.RDS")

saveRDS(prov_dt, "data/shapes/simadmin.RDS")

saveRDS(spain_dt, "data/shapes/spainshape.RDS")

saveRDS(incomedata, "data/incomedata.RDS")















