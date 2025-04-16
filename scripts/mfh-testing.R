################################################################################
########## IMPLEMENT THE ISABEL MSAE PACKAGE FOR SMALL AREA ESTIMATION #########
################################################################################

### load packages

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(msae, tidyverse, data.table, sae, car)


### include a poverty line and calculate some indicators
data("incomedata")

incomedata <- 
  incomedata %>%
  mutate(povline = 6557.143) %>%
  mutate(y = (povline - income) / povline) %>%
  mutate(poverty = ifelse(y > 0, TRUE, FALSE),
         y0 = ifelse(y > 0, 1, 0),
         y1 = ifelse(y > 0, y, 0),
         y2 = ifelse(y > 0, y^2, 0))

### let create some additional variables
incomedata <- 
  incomedata %>%
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

### estimated domain size
sampsize <- aggregate(incomedata$weight, list(incomedata$prov), sum)[, 2]

### estimate P0, P1 and P2
prov_est <- 
  incomedata %>%
  group_by(prov) %>%
  summarise(p0_prov = sum(weight * y0),
            p1_prov = sum(weight * y1),
            p2_prov = sum(weight * y2)) %>%
  mutate(p0_prov = p0_prov / sampsize,
         p1_prov = p1_prov / sampsize,
         p2_prov = p2_prov / sampsize)

incomedata <- 
  incomedata %>% left_join(prov_est, by = "prov")


prov_variance <- 
  incomedata %>%
  mutate(v11 = ifelse(poverty, (weight * (weight - 1)) * (y0 - p0_prov)*(y0 - p0_prov), 0),
         v12 = ifelse(poverty, (weight * (weight - 1)) * (y0 - p0_prov)*(y1 - p1_prov), 0),
         v13 = ifelse(poverty, (weight * (weight - 1)) * (y0 - p0_prov)*(y2 - p2_prov), 0),
         v22 = ifelse(poverty, (weight * (weight - 1)) * (y1 - p1_prov)*(y1 - p1_prov), 0),
         v23 = ifelse(poverty, (weight * (weight - 1)) * (y1 - p1_prov)*(y2 - p2_prov), 0),
         v33 = ifelse(poverty, (weight * (weight - 1)) * (y2 - p2_prov)*(y2 - p2_prov), 0)) %>%
  group_by(prov) %>%
  summarise_at(c("v11", "v12", "v13", "v22", "v23", "v33"), sum) %>%
  mutate_at(c("v11", "v12", "v13", "v22", "v23", "v33"),
            function(x){ 
              
              x / sampsize^2
              
              })


#### quickly carry out variable selection for each outcome variable

candidate_vars <- colnames(incomedata)[!colnames(incomedata) %in% 
                                         c("provlab", "prov", "income",
                                           "weight", "povline", "y",
                                           "poverty", "y0", "y1", "y2",
                                           "p0_prov", "p1_prov", "p2_prov",
                                           "ac", "nat", "educ", "labor")]

prov_dt <- 
  incomedata %>%
  as.data.table() %>%
  .[, lapply(.SD, weighted.mean, w = weight, na.rm = TRUE),
    .SDcols = candidate_vars,
    by = c("provlab", "prov")]
  

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



prov_dt <- merge(prov_dt, prov_est, by = "prov")

### select the set of candidate variables by dropping the 
### very highly correlated variables

outcome_list <- c("p0_prov", "p1_prov", "p2_prov")


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
  
vardir <- c("v11", "v12", "v13", "v22", "v23", "v33")


### now lets estimate all 4 models
prov_dt <- merge(prov_dt, prov_variance, by = "prov")

model0_obj <- eblupUFH(mfh_formula, vardir = vardir, data = prov_dt)
model1_obj <- eblupMFH1(mfh_formula, vardir = vardir, data = prov_dt)
model2_obj <- eblupMFH2(mfh_formula, vardir = vardir, data = prov_dt)
# model3_obj <- eblupMFH3(mfh_formula, vardir = vardir, data = prov_dt, PRECISION = 1e-4)


#### now we need to check all the assumptions















