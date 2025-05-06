################################################################################
########### ESTIMATING THE MFH SPATIO TEMPORAL MODEL WITH SPAIN DATA ###########
################################################################################

### load the libraries

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(sf, data.table, tidyverse, car, msae, sae)


### read in the datasets

income_dt <- readRDS("data/incomedata.RDS")

shp_dt <- readRDS("data/shapes/spainshape.RDS")


##### let us start by developing the typical MFH model 

#### step 1: compute the direct estimates for the poverty rates

income_dt <- 
  income_dt %>%
  mutate(y2012 = (povline2012 - income2012) / povline2012,
         y2013 = (povline2013 - income2013) / povline2013,
         y2014 = (povline2014 - income2014) / povline2014) %>%
  mutate(poverty2012 = ifelse(y2012 > 0, TRUE, FALSE),
         poverty2013 = ifelse(y2013 > 0, TRUE, FALSE),
         poverty2014 = ifelse(y2014 > 0, TRUE, FALSE),
         y2012 = ifelse(y2012 > 0, 1, 0),
         y2013 = ifelse(y2013 > 0, 1, 0),
         y2014 = ifelse(y2014 > 0, 1, 0))


### let create some additional variables
income_dt <- 
  income_dt %>%
  mutate(across(c(gen, nat), ~ case_when(
    .x == 1 ~ 0,
    .x == 2 ~ 1,
    TRUE ~ NA_real_
  )))

### create some quick variable interactions
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

### estimated domain size
sampsize <- aggregate(income_dt$weight, list(income_dt$prov), sum)[, 2]


#' A function to compute variance-covariance matrix from unit level data 
#' 
#' This function will compute the variance covariance matrix from the unit level data. 
#' E.g. household survey with income/expenditure variable to be used for poverty mapping
#' 
#' 
#' @param dt a `data.frame` object i.e. the household survey 
#' @param y a `character` vector, a list of outcome variables of interest. 
#' @param weights a `character`, weight variable name, default is NULL
#' @param povline a `numeric`, the poverty line. Could be one value applied to each 
#' element of `y` or equal in length to `y`
#' @param domain `character`, the target area variable name
#' @export
#' 

compute_covar_matrix <- function(dt, y, weights, povline, domain){
  
  dt <- 
    dt %>%
    dplyr::select(y, weights, domain) %>%
    as.data.frame()
  
  ifelseworker <- function(x, y){
    
    zz <- ifelse(x < y, 1, 0)
    
    return(zz)
    
  }
  
  if (length(povline) == 1 & length(y) > 1){
    
    povline <- rep(povline, length(y))
    
  }
  
  ind_dt <- 
    purrr::map2(.x = dt[, y],
                .y = povline,
                .f = ~ ifelseworker(.x, .y)) %>%
    as.data.frame() %>%
    cbind(dt %>% dplyr::select(domain, weights)) %>%
    as.data.table()
  
  ind_dt[, paste0("avg", y) :=
         lapply(.SD, function(x) weighted.mean(x, w = get(weights))), 
         .SDcols = y, 
         by = domain]
    
  ### compute the variance-covariance matrices
  ##### create the appropriate list of pairs
  pair_list <- c(lapply(y, rep, 2),
                 combn(y, 2, simplify = FALSE))
  
  ind_dt <- 
  ind_dt %>%
    rename(weights = weights)
  
  
  popsize <- aggregate(ind_dt$weights,
                       list(ind_dt[[domain]]),
                       sum)[,2]


  ### computation time!
  comp_covar <- function(vars){
    
    ind_dt <- 
      ind_dt %>%
      mutate(ind1 = !!sym(vars[[1]]),
             ind2 = !!sym(vars[[2]]),
             avg1 = !!sym(paste0("avg", vars[[1]])),
             avg2 = !!sym(paste0("avg", vars[[2]])))

    prov_dt <- 
      ind_dt %>%
      mutate(v = ifelse(ind1 == 1, 
                        (weights * (ind1 - avg1) * (ind2 - avg2)), 
                        0)) %>%
      group_by(!!sym(domain)) %>%
      summarize(v = sum(v, na.rm = TRUE)) %>%
      mutate(v = v / popsize)
    
    return(prov_dt)

  }
  
  var_dt <- 
    lapply(X = pair_list,
           comp_covar)
  
  ### construct names
  pair_list <- lapply(pair_list,
                      function(x){
                        
                        zz <- paste0("v_", paste(x, collapse = ""))
                        
                        return(zz)
                        
                      })
  var_dt <- 
    mapply(FUN = function(dt, y){
      
      dt <- dt %>% rename(!!sym(y) := v)
      
      return(dt)
      
    }, SIMPLIFY = FALSE,
    dt = var_dt,
    y = pair_list)
  
  var_dt <- Reduce(f = merge,
                   x = var_dt)
   
  return(var_dt)
  
  
}


### include province poverty rates 
povcols <- c("y2012", "y2013", "y2014")

income_dt <- 
  income_dt %>%
  as.data.table() %>%
  .[, paste0(povcols, "prov") := lapply(.SD, weighted.mean, w = weight),
    .SDcols = povcols,
    by = prov]


var_dt <- cov.wt(x = income_dt[, c("y2012", "y2013", "y2014")],
                 wt = income_dt$weight)


# ### now lets include the covariance matrix as well
# var_dt <- compute_covar_matrix(dt = income_dt,
#                                y = c("income2012", "income2013", "income2014"),
#                                weights = "weight",
#                                povline = c(unique(income_dt$povline2012),
#                                            unique(income_dt$povline2013),
#                                            unique(income_dt$povline2014)),
#                                domain = "prov")
# 
# income_dt[, sampsize := sum(weight, na.rm = TRUE), by = prov]

var_dt <- as.numeric(c(diag(var_dt$cov), var_dt$cov[lower.tri(var_dt$cov, diag = FALSE)]))

names(var_dt) <- c("v1", "v2", "v3", "v12", "v13", "v23")

### create the province level dataset

candidate_vars <- colnames(income_dt)[!colnames(income_dt) %in% 
                                         c("provlab", "prov", "income",
                                           "weight", "povline", "y",
                                           "poverty", "y0", "y1", "y2",
                                           "p0_prov", "p1_prov", "p2_prov",
                                           "ac", "nat", "educ", "labor",
                                           "age")]

candidate_vars <- candidate_vars[!grepl("sampsize|prov|poverty|income|povline|^v[0-9]|^y[0-9]", 
                                        candidate_vars)]


prov_dt <- 
  income_dt[, lapply(.SD, weighted.mean, w = weight), 
                     .SDcols = c(candidate_vars, "y2012", "y2013", "y2014"),
                     by = "prov"] %>%
  mutate(!!!as.list(var_dt))

  
#### lets do some variable selection
### a quick function to perform stepwise variable selection both ways
#' A function to perform stepwise variable selection with AIC selection criteria
#' 
#' @param dt data.frame, dataset containing the set of outcome and independent variables
#' @param xvars character vector, the set of x variables
#' @param y chr, the name of the y variable
#' @param vif_threshold integer, a variance inflation factor threshold
#' 
#' @import data.table
#' @importFrom cars vif
#' @importFrom MASS stepAIC

stepAIC_wrapper <- function(dt, xvars, y) {
  
  dt <- as.data.table(dt)
  
  # Drop columns that are entirely NA
  dt <- dt[, which(unlist(lapply(dt, function(x) !all(is.na(x))))), with = FALSE]
  
  xvars <- xvars[xvars %in% colnames(dt)]
  
  # Keep only complete cases
  dt <- na.omit(dt[, c(y, xvars), with = FALSE])
  
  # Drop multicollinear variables
  model_formula <- as.formula(paste(y, " ~ ", paste(xvars, collapse = " + ")))
  lm_model <- lm(model_formula, data = dt)
  
  # Check for aliased coefficients
  aliased <- is.na(coef(lm_model))
  if (any(aliased)) {
    xvars <- names(aliased)[!aliased & names(aliased) != "(Intercept)"]
  }
  
  # Re-fit after removing aliased
  model_formula <- as.formula(paste(y, "~", paste(xvars, collapse = " + ")))
  full_model <- lm(model_formula, data = dt)
  
  # Apply stepwise selection
  stepwise_model <- stepAIC(full_model, direction = "both", trace = 0)
  
  return(stepwise_model)
  
}




### select the set of candidate variables by dropping the 
### very highly correlated variables

outcome_list <- c("y2012", "y2013", "y2014")

candidate_vars <- candidate_vars[!candidate_vars %in% c("yos", "abc", "poi",
                                                        "sam", "ntl", "ips")]

stepaicmodel_list <- 
  mapply(x = outcome_list,
         # y = selvar_list,
         FUN = function(x){
           
           y <- stepAIC_wrapper(dt = prov_dt,
                                xvars = candidate_vars,
                                y = x)
           
           coef_list <- y[["coefficients"]]
           
           coef_list <- names(coef_list)[!names(coef_list) %in% "(Intercept)"]
           
           return(coef_list)
           
         },
         SIMPLIFY = FALSE)



### create the formulas

mfh_formula <-
  mapply(x = outcome_list,
         y = stepaicmodel_list,
         FUN = function(x, y){
           
           formula <- as.formula(paste0(x, " ~ ", paste(y, collapse = " + "))) 
           
           return(formula)
           
         }, SIMPLIFY = FALSE)



### now lets estimate all 4 models


model0_obj <- eblupUFH(mfh_formula, vardir = names(var_dt), data = prov_dt)
model1_obj <- eblupMFH1(mfh_formula, vardir = names(var_dt), data = prov_dt)
model2_obj <- eblupMFH2(mfh_formula, vardir = names(var_dt), data = prov_dt, MAXITER = 10000)
model3_obj <- eblupMFH3(mfh_formula, vardir = names(var_dt), data = prov_dt, MAXITER = 10000)


### lets check the normality of the random effects






















