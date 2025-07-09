#' A function to perform stepwise variable selection with AIC selection criteria
#' 
#' @param dt `data.frame`, dataset containing the set of outcome and independent variables
#' @param xvars `character`, the character vector of variable names
#' @param y `character`, the name of the y variable
#' @param cor_threshold `numeric`, value between 0 and 1 for the correlation threshold beyond 
#' which variable should be dropped
#' @param k `numeric` the multiple of the number of degrees of freedom used for the penalty. 
#' Only k = 2 gives the genuine AIC: k = log(n) is sometimes referred to as BIC or SBC.
#' 
#' @import data.table
#' @importFrom cars vif
#' @importFrom MASS stepAIC

step_wrapper <- function(dt, 
                         xvars, 
                         y, 
                         cor_thresh = 0.95, 
                         k) {
  
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
  stepwise_model <- stepAIC(full_model, 
                            direction = "both", 
                            trace = 0,
                            k = k)
  
  return(stepwise_model)
  
}
