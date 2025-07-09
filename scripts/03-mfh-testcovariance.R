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


### direct estimation within the fay herriot context
data(sizeprov)

## quickly compute sample size for each province
sampsize_dt <- 
  income_dt |>
  group_by(prov) |>
  summarize(N = n())

## the poverty line for each is already included within the data. 
## Lets compute the direct estimates for each year of data and create a list of data.frames (equal in length to the number of years) 
## containing the direct estimate, the standard errors and the coefficient of variation. 

direct_list <- 
  mapply(FUN = function(y, threshold){
    
    z <- emdi::direct(y = y,
                      smp_data = income_dt %>% as.data.table(),
                      smp_domains = "prov",
                      weights = "weight",
                      threshold = unique(income_dt[[threshold]]),
                      var = TRUE)
    
    z <- 
      z$ind |>
      dplyr::select(Domain, Head_Count) |>
      rename(Direct = "Head_Count") |>
      merge(z$MSE |>
              dplyr::select(Domain, Head_Count) |>
              rename(MSE = "Head_Count"),
            by = "Domain") |>
      mutate(SD = sqrt(MSE)) |>
      mutate(CV = SD / Direct) |>
      merge(sampsize_dt |> 
              mutate(prov = as.factor(prov)), 
            by.x = "Domain", 
            by.y = "prov") |>
      mutate(var_SRS = Direct * (1 - Direct) / N) |>
      mutate(deff = MSE / var_SRS) |>
      mutate(n_eff = N / deff)
    
    
    return(z)
    
  }, SIMPLIFY = FALSE,
  y = c("income2012", "income2013", "income2014"),
  threshold = c("povline2012", "povline2013", "povline2014"))



## quickly creating the poverty indicator
income_dt <- 
  income_dt |>
  mutate(poor2012 = ifelse(income2012 < povline2012, 1, 0),
         poor2013 = ifelse(income2013 < povline2013, 1, 0),
         poor2014 = ifelse(income2014 < povline2014, 1, 0))


### a simple worker function for computing the variance covariance matrix of the sampling error of the direct estimate

### running the compute_vcov function to prepare the variance covariance matrix 
var_dt <- 
  compute_vcov(dt = income_dt, 
               domain = "prov",
               ids = 1, 
               weights = "weight", 
               yvars = c("poor2012", "poor2013", "poor2014"))

### sample size
var_dt <- 
  var_dt |> 
  merge(sampsize_dt, by.x = "domain", by.y = "prov")



vardir <- grep("^v_", names(var_dt), value = TRUE)

var_dt <-
  lapply(X = vardir,
         FUN = function(x){
           
           z <- varsmoothie_king(domain = var_dt[["domain"]],
                                 direct_var = var_dt[[x]],
                                 sampsize = var_dt[["N"]]) |>
             as.data.table() |>
             setnames(old = "var_smooth", new = paste0("vs", x)) |>
             as_tibble()
           
           
           return(z)
           
         }) %>%
  Reduce(f = "merge",
         x = .) %>%
  merge(x = var_dt,
        y = .,
        by.x = "domain",
        by.y = "Domain") |>
  as_tibble()



### variable selection
### the set of outcome variables
candidate_vars <- colnames(prov_dt)[!colnames(prov_dt) %in% c("prov", "provlab")]


### creating the province level data for variable selection and poverty mapping
prov_dt <- 
  map2(.x = list("poor2012", "poor2013", "poor2014"), ### first we create a direct estimate dataset from direct_list
       .y = direct_list,
       ~ .y |> 
         rename(!!sym(.x) := Direct) |>
         dplyr::select(Domain, all_of(.x))) |>
  Reduce(f = "merge") |> ### merge the selected variables from each dataframe of the list
  merge(prov_dt, ### merging the prepared simulated administrative data variables
        by.x = "Domain", by.y = "prov", all = TRUE) |>
  merge(var_dt,
        by.x = "Domain", by.y = "domain", all = TRUE)

outcome_list <- list("poor2012", "poor2013", "poor2014")

fh_step <-
  lapply(X = outcome_list,
         FUN = function(x){
           
           model_obj <-
             step_wrapper(dt = prov_dt,
                          xvars = candidate_vars,
                          y = x,
                          cor_thresh = 0.7,
                          k = log(nrow(prov_dt))) ### using log(n) to force BIC selection
           
           xx <- names(model_obj$coefficients)[!grepl("(Intercept)",
                                                      names(model_obj$coefficients))]
           return(xx)
           
         })

# Then generate the list of formulas
mfh_formula <- 
  mapply(x = fh_step,
         y = outcome_list,
         function(x, y) {
           
           as.formula(paste(y, " ~ ", paste(x, collapse = " + "))) 
           
         }, SIMPLIFY = FALSE)

names(mfh_formula) <- outcome_list


lmcheck_obj <- 
  lapply(X = mfh_formula,
         FUN = lm,
         data = prov_dt)


lapply(X = lmcheck_obj,
       FUN = summary)

### variance-covariance matrix columns
varcols <- colnames(var_dt)[grepl(pattern = "^v_", x = colnames(var_dt))] 

## replace the variances-covariances that are zero with their smoothed counterparts
prov_dt <- 
  prov_dt |>
  mutate(across(
    starts_with("v_"),
    ~ if_else(abs(.x) <= 1e-4, get(paste0("vsv", str_remove(cur_column(), "^v"))), .x),
    .names = "{.col}"
  ))



#### now we estimate all 4 models including the multivariate fay herriot models

model0_obj <- eblupUFH(mfh_formula, vardir = varcols, data = prov_dt)
model1_obj <- eblupMFH1(mfh_formula, vardir = varcols, data = prov_dt, MAXITER = 1e10, PRECISION = 1e-2)
model2_obj <- eblupMFH2(mfh_formula, vardir = varcols, data = prov_dt, MAXITER = 1e10, PRECISION = 1e-2)
# model3_obj <- eblupMFH3(mfh_formula, vardir = varcols, data = prov_dt, MAXITER = 1e10, PRECISION = 1e-2)



simulation <- function(prov_dt,
                       MAXITER = 1e10,
                       PRECISION = 1e-2){
  
    
  #### now we estimate all 4 models including the multivariate fay herriot models
  model0_obj <- eblupUFH(mfh_formula, 
                         vardir = varcols, 
                         data = prov_dt)
  
  model1_obj <- eblupMFH1(mfh_formula, 
                          vardir = varcols, 
                          data = prov_dt, 
                          MAXITER = MAXITER, 
                          PRECISION = PRECISION)
  
  model2_obj <- eblupMFH2(mfh_formula, 
                          vardir = varcols, 
                          data = prov_dt, 
                          MAXITER = MAXITER, 
                          PRECISION = PRECISION)
  
  
  ### lets create a dataframe with the errors and estimated poverty rates
  
  eblup_dt <- model2_obj$eblup
  mse_dt <- model2_obj$MSE
  
  colnames(eblup_dt) <- paste0("eblup_", colnames(eblup_dt))
  colnames(mse_dt) <- paste0("mse_", colnames(mse_dt))
  
  reset_dt <- bind_cols(eblup_dt, mse_dt) |> as_tibble()
  
  ### lets perform the reset test on the pairs of variables as appropriate i.e. poverty = B0 + B1*MSE
  
  test_list <- 
    lapply(outcome_list, 
           function(x) {
      
              # Subset only the columns for year x
              dt <- 
                reset_dt |>
                dplyr::select(matches(paste0(x, "$")))  # Select columns ending with current year
              
              yvar <- colnames(dt)[grepl("^eblup_", colnames(dt))]
              xvar <- colnames(dt)[grepl("^mse_", colnames(dt))]
              
              # Create formula with squared and cubed terms using I()
              form <- as.formula(paste0(
                xvar, " ~ ", 
                yvar, " + I(", yvar, "^2)"
              ))
              
              model_obj <- lm(form, data = dt)
              
              return(model_obj)
      
         })
  
  reset_dt2 <- 
    reset_dt |> 
    mutate(id = row_number()) |> 
    pivot_longer(c(contains("eblup"), contains("mse"))) |> 
    mutate(
      year = str_extract(name, "\\d{4}") |> as.numeric(),
      var =  str_extract(name, "eblup|mse")
    ) |> 
    dplyr::select(-name) |> 
    pivot_wider(names_from = var, values_from = value )
  
  reset_dt2 |> 
    ggplot() + 
    aes(x = eblup, y = mse) + 
    geom_point() + 
    facet_wrap(. ~ year) + 
    geom_smooth() + 
    theme_bw()
  
  
  
  
}

































