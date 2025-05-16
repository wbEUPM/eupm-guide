################################################################################
########### ESTIMATING THE MFH SPATIO TEMPORAL MODEL WITH SPAIN DATA ###########
################################################################################

### load the libraries

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(sf, data.table, tidyverse, car, msae, sae, survey, spdep)


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


# var_dt <- cov.wt(x = income_dt[, c("y2012", "y2013", "y2014")],
#                  wt = income_dt$weight)
# 
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

design_obj <- svydesign(ids = ~1, weights = ~weight, data = income_dt)

var_dt <- svymean(~y2012 + y2013 + y2014, design_obj)

var_dt <- vcov(var_dt)

var_dt <- as.numeric(c(diag(var_dt), var_dt[lower.tri(var_dt, diag = FALSE)]))

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

stepAIC_wrapper <- function(dt, xvars, y, cor_thresh = 0.95) {
  
  library(data.table)
  library(MASS)        # For stepAIC
  library(caret)       # For findLinearCombos
  
  dt <- as.data.table(dt)
  
  # Drop columns that are entirely NA
  dt <- dt[, which(unlist(lapply(dt, function(x) !all(is.na(x))))), with = FALSE]
  
  xvars <- xvars[xvars %in% colnames(dt)]
  
  # Keep only complete cases
  dt <- na.omit(dt[, c(y, xvars), with = FALSE])
  
  # Step 1: Remove aliased (perfectly collinear) variables
  model_formula <- as.formula(paste(y, "~", paste(xvars, collapse = " + ")))
  lm_model <- lm(model_formula, data = dt)
  aliased <- is.na(coef(lm_model))
  if (any(aliased)) {
    xvars <- names(aliased)[!aliased & names(aliased) != "(Intercept)"]
  }
  
  # Step 2: Remove near-linear combinations
  xmat <- as.matrix(dt[, ..xvars])
  combo_check <- tryCatch(findLinearCombos(xmat), error = function(e) NULL)
  if (!is.null(combo_check) && length(combo_check$remove) > 0) {
    xvars <- xvars[-combo_check$remove]
    xmat <- as.matrix(dt[, ..xvars])
  }
  
  # Step 3: Drop highly correlated variables
  cor_mat <- abs(cor(xmat))
  diag(cor_mat) <- 0
  while (any(cor_mat > cor_thresh, na.rm = TRUE)) {
    cor_pairs <- which(cor_mat == max(cor_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
    var1 <- colnames(cor_mat)[cor_pairs[1]]
    var2 <- colnames(cor_mat)[cor_pairs[2]]
    # Drop the variable with higher mean correlation
    drop_var <- if (mean(cor_mat[var1, ]) > mean(cor_mat[var2, ])) var1 else var2
    xvars <- setdiff(xvars, drop_var)
    xmat <- as.matrix(dt[, ..xvars])
    cor_mat <- abs(cor(xmat))
    diag(cor_mat) <- 0
  }
  
  # Step 4: Warn if still ill-conditioned
  cond_number <- kappa(xmat, exact = TRUE)
  if (cond_number > 1e10) {
    warning("Design matrix is ill-conditioned (condition number > 1e10). Consider reviewing variable selection.")
  }
  
  # Final model fit
  model_formula <- as.formula(paste(y, "~", paste(xvars, collapse = " + ")))
  full_model <- lm(model_formula, data = dt)
  
  # Stepwise selection
  stepwise_model <- stepAIC(full_model, direction = "both", trace = 0)
  
  return(stepwise_model)
  
}



### select the set of candidate variables by dropping the 
### very highly correlated variables

outcome_list <- c("y2012", "y2013", "y2014")

# candidate_vars <- candidate_vars[!candidate_vars %in% c("yos", "abc", "poi",
#                                                         "sam", "ntl", "ips")]

stepaicmodel_list <- 
  mapply(x = outcome_list,
         # y = selvar_list,
         FUN = function(x){
           
           y <- stepAIC_wrapper(dt = prov_dt,
                                xvars = candidate_vars,
                                y = x,
                                cor_thresh = 0.8)
           
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
model3_obj <- eblupMFH3(mfh_formula, vardir = names(var_dt), data = prov_dt, MAXITER = 1e7, PRECISION = 1e-3)

saveRDS(model3_obj, "data/modelmfh3.RDS")


### lets check the normality of the random effects
lmcheck_obj <- 
lapply(X = mfh_formula,
       FUN = lm,
       data = prov_dt)


lapply(X = lmcheck_obj,
       FUN = summary)


### evaluating the normality assumption
resid_dt <- prov_dt[,c("y2012", "y2013", "y2014")] - as.data.table(model3_obj$eblup)

### perform the shapiro test

shapiro_obj <- apply(resid_dt, 2, shapiro.test)

summary_dt <- 
  data.frame(Time = names(shapiro_obj),
             W = lapply(X = shapiro_obj,
                        FUN = function(x){
                          
                          return(x$statistic[[1]])
                          
                        }) %>%
               as.numeric(),
             p_value = lapply(X = shapiro_obj,
                              FUN = function(x){
                                
                                return(x$p.value)
                                
                              }) %>%
               as.numeric())

### plot the results
summary_dt <- 
  summary_dt %>%
  mutate(label = paste0("W = ", round(W, 3), "\n", "p = ", signif(p_value, 3)))

resid_dt %>%
  pivot_longer(cols = everything(), 
               names_to = "Time", 
               values_to = "Residual") %>%
  ggplot(aes(x = Residual)) + 
  geom_histogram(bins = 10, fill = "steelblue", color = "white") + 
  geom_text(data = summary_dt, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~Time, scales = "free") + 
  theme_minimal() + 
  labs(title = "Residual Histograms by Time Period")


### here's how to create qqplots
resid_dt %>%
  pivot_longer(cols = everything(),
               names_to = "Time",
               values_to = "Residual") %>%
  ggplot(aes(sample = Residual)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~Time, scales = "free") +
  theme_minimal() +
  labs(title = "QQ Plots of Residuals by Time Period")



#### For the random effects
raneff_dt <- as.data.frame(model3_obj$randomEffect)

### lets run the shapiro wilks tests again
shapiro_obj <- apply(raneff_dt, 2, shapiro.test)


summary_dt <- 
  data.frame(Time = names(shapiro_obj),
             W = lapply(X = shapiro_obj,
                        FUN = function(x){
                          
                          return(x$statistic[[1]])
                          
                        }) %>%
               as.numeric(),
             p_value = lapply(X = shapiro_obj,
                              FUN = function(x){
                                
                                return(x$p.value)
                                
                              }) %>%
               as.numeric())

### plot the results
summary_dt <- 
  summary_dt %>%
  mutate(label = paste0("W = ", round(W, 3), "\n", "p = ", signif(p_value, 3)))

raneff_dt %>%
  pivot_longer(cols = everything(), 
               names_to = "Time", 
               values_to = "RandEff") %>%
  ggplot(aes(x = RandEff)) + 
  geom_histogram(bins = 10, fill = "darkorange", color = "white") + 
  geom_text(data = summary_dt, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~Time, scales = "free") + 
  theme_minimal() + 
  labs(title = "Randon Effects Histograms by Time Period")

### lets do spatio-temporal (we used https://paezha.github.io/spatial-analysis-r/area-data-ii.html)
shp_dt <- 
  shp_dt %>% 
  sf::st_as_sf(crs = 4326,
               agr = "constant")

### create a neighborhood list based on queen continguity
shp_nb <- poly2nb(pl = shp_dt, 
                  queen = TRUE, 
                  row.names = shp_dt$provlab %>% 
                    st_drop_geometry() %>%
                    as.character())

### convert to a spatial weights matrix
lw <- nb2mat(shp_nb, style = "B", zero.policy = TRUE)

colnames(lw) <- 
  shp_dt$provlab %>% 
  st_drop_geometry() %>% 
  as.character()

# Find provinces with all-zero rows and/or columns
islands <- union(
  rownames(lw)[rowSums(lw) == 0],
  colnames(lw)[colSums(lw) == 0]
)

# Create full pairings between island provinces and all others
island_pairs <- 
  expand.grid(Var1 = islands, Var2 = colnames(lw)) %>%
  bind_rows(expand.grid(Var1 = rownames(lw), Var2 = islands)) %>%
  distinct()


dist_dt <- 
  lw %>%
  as.matrix() %>%
  as.table() %>%
  as.data.frame() %>%
  filter(Freq != 0) %>%
  dplyr::select(Var1, Var2) %>%
  bind_rows(island_pairs)



dist_list <- 
    lapply(1:nrow(dist_dt), 
           function(i) {
                        
      pair <- as.character(unlist(dist_dt[i,]))
      
      dist <- shp_dt %>%
        dplyr::filter(provlab %in% pair) %>%
        st_centroid() %>%
        { st_distance(.[1, ], .[2, ]) } %>%
        as.numeric()
      
      return(dist)
    
  })


dist_dt$dist <- unlist(dist_list)

### reinclude this within the lw matrix and then ensure each row sums up to 1

lw_dt <- 
  lw %>%
  as.table() %>%
  as.data.frame() %>%
  dplyr::select(Var1, Var2) %>%
  merge(dist_dt, all.x = TRUE) %>%
  mutate(dist = ifelse(is.na(dist), 0, dist)) %>%
  pivot_wider(names_from = Var2, values_from = dist) %>%
  column_to_rownames(var = "Var1") %>%
  apply(X = ., 
        MARGIN = 1, 
        FUN = function(x){
          
          x <- as.numeric(x)
          
          z <- x / sum(x)
          
          return(z)
          
        }) %>%
  as.matrix(dimnames = list(dist_dt$Var1, dist_dt$Var2))



#### reprepare the structure of the data to fit
long_dt <- 
  prov_dt[, c("y2012", "y2013", "y2014", 
              candidate_vars, "prov"), 
          with = F] %>%
  pivot_longer(cols = starts_with("y"),
               names_to = "year",
               names_prefix = "y",
               values_to = "value") %>%
  mutate(year = as.integer(year))


### compute the variance 
var_dt <- 
  mapply(y = c("income2012", "income2013", "income2014"),
         threshold = c("povline2012", "povline2013", "povline2014"),
         FUN = function(y, threshold){
           
           
           z <- emdi::direct(y = y,
                             smp_data = income_dt,
                             smp_domains = "prov",
                             weights = "weight",
                             threshold = unique(income_dt[[threshold]]),
                             var = TRUE)
           
           z <- 
             z$MSE[, c("Domain", "Head_Count")] %>%
             rename(prov = "Domain",
                    var = "Head_Count") %>%
             mutate(year = threshold)
           
           
           
           return(z)
           
         }, 
         SIMPLIFY = FALSE)

var_dt <- Reduce(f = "rbind",
                 x = var_dt)

var_dt <- 
  var_dt %>% 
  mutate(year = as.integer(substr(year, 
                                  nchar(year) - 3, 
                                  nchar(year))))

long_dt <- merge(long_dt, var_dt, by = c("prov", "year"))

rownames(lw_dt) <- colnames(lw_dt)

selvars_list <- 
stepAIC_wrapper(dt = long_dt,
                xvars = candidate_vars,
                y = "value")

selvars_list <- names(selvars_list$coefficients)[!names(selvars_list$coefficients) %in% "(Intercept)"]


stfh_formula <- as.formula(paste0("value", " ~ ", paste(selvars_list, collapse = " + ")))


stfh_obj <- pbmseSTFH(formula = stfh_formula,
                      D = nrow(lw_dt),
                      T = length(unique(long_dt$year)),
                      vardir = var,
                      proxmat = lw_dt,
                      data = long_dt)


stfh_cv <- 100 * stfh_obj$mse / stfh_obj$est$eblup

direct_cv <- 100 * sqrt(long_dt$var) / long_dt$value

results_dt <- data.frame(area = long_dt$prov,
                         time = long_dt$year,
                         direct = long_dt$value,
                         eblup_stfh = stfh_obj$est$eblup,
                         direct_cv = sqrt(long_dt$var) / long_dt$value,
                         stfh_cv = sqrt(stfh_obj$mse) / stfh_obj$est$eblup)


### comparing the direct vs the eblup estimates
results_dt %>%
  pivot_longer(cols = c("direct", "eblup_stfh"),
               names_to = "method",
               values_to = "estimate") %>%
  ggplot(aes(x = factor(area), 
             y = estimate, 
             color = method, 
             group = method)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~time, ncol = 1) +
  labs(
    x = "Area",
    y = "Estimate",
    title = "Direct vs EBLUP Estimates by Area and Time",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


### comparing the coefficient of variation (reduction in variability)
results_dt %>%
  pivot_longer(cols = c("direct_cv", "stfh_cv"),
               names_to = "method",
               values_to = "estimate") %>%
  ggplot(aes(x = factor(area), 
             y = estimate, 
             color = method, 
             group = method)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~time, ncol = 1) +
  labs(
    x = "Area",
    y = "CV",
    title = "Direct vs EBLUP CVs by Area and Time",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


### comparing the stfh, univariate FH and the direct estimation
ufh_list <- 
long_dt %>%
  stepAIC_wrapper(dt = .,
                  xvars = candidate_vars,
                  y = "value") %>%
  .[["coefficients"]] %>%
  names() %>%
  setdiff("(Intercept)") 

ufh_formula <- as.formula(paste0("value", " ~ ", paste(ufh_list, collapse = " + ")))

ufh_time1 <- mseFH(ufh_formula, vardir = var, data = long_dt[long_dt$year == 2012,])
ufh_time2 <- mseFH(ufh_formula, vardir = var, data = long_dt[long_dt$year == 2013,])
ufh_time3 <- mseFH(ufh_formula, vardir = var, data = long_dt[long_dt$year == 2014,])

### create the table of fh estimates and cvs
fh_dt <- rbind(data.frame(area = long_dt$prov[long_dt$year == 2012], 
                          time = long_dt$year[long_dt$year == 2012], 
                          eblup_fh = ufh_time2$est$eblup, 
                          fh_cv = sqrt(ufh_time2$mse) / ufh_time2$est$eblup),
               data.frame(area = long_dt$prov[long_dt$year == 2013], 
                          time = long_dt$year[long_dt$year == 2013], 
                          eblup_fh = ufh_time2$est$eblup, 
                          fh_cv = sqrt(ufh_time2$mse) / ufh_time2$est$eblup),
               data.frame(area = long_dt$prov[long_dt$year == 2014], 
                          time = long_dt$year[long_dt$year == 2014], 
                          eblup_fh = ufh_time3$est$eblup, 
                          fh_cv = sqrt(ufh_time3$mse) / ufh_time3$est$eblup))

plot_dt <- 
  rbind(results_dt %>%
          pivot_longer(cols = c("direct", "eblup_stfh"),
                       names_to = "method",
                       values_to = "estimate") %>%
          dplyr::select(area, time, estimate, method),
        fh_dt %>% 
          pivot_longer(cols = c("eblup_fh"),
                       names_to = "method",
                       values_to = "estimate") %>%
          dplyr::select(area, time, estimate, method))


plot_dt %>%
  ggplot(aes(x = as.factor(time), 
             y = estimate, 
             color = method, 
             group = method)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~area, ncol = 3, scales = "free_y") +
  labs(
    x = "Time",
    y = "Estimates",
    title = "Comparing Stability of Spatio-Temporal FH vs FH vs Direct Estimates",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



### now lets show that stfh estimates are more stable than FH and direct over time
plot_dt <- 
  rbind(results_dt %>%
          pivot_longer(cols = c("direct_cv", "stfh_cv"),
                       names_to = "method",
                       values_to = "estimate") %>%
          dplyr::select(area, time, estimate, method),
        fh_dt %>% 
          pivot_longer(cols = c("fh_cv"),
                       names_to = "method",
                       values_to = "estimate") %>%
          dplyr::select(area, time, estimate, method))

plot_dt[plot_dt$area <= 20,] %>%
  ggplot(aes(x = as.factor(time), 
             y = estimate, 
             color = method, 
             group = method)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~area, ncol = 3, scales = "free_y") +
  labs(
    x = "Time",
    y = "CV",
    title = "Comparing Stability of Spatio-Temporal FH vs FH vs Direct Estimates",
    color = "Method"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


list.files(path = "scripts",
           pattern = "^Function_.*\\.R$",
           full.names = TRUE) %>%
  lapply(X = .,
         FUN = source)


















































