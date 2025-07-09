#' A function to perform stepwise variable selection based on selection criteria
#' 
#' @param dt data.frame, dataset containing the set of outcome and independent variables
#' @param xvars character vector, the set of x variables
#' @param y character, the name of the y variable
#' @param cor_thresh double, a correlation threshold between 0 and 1
#' @param criteria character, criteria that can be chosen are "AIC", "AICc", "AICb1", "AICb2", "BIC", "KIC", "KICc", "KICb1", or "KICb2". Defaults to "AIC". If transformation is set to "arcsin", only "AIC" and "BIC" can be chosen.
#' @param vardir character, name of the variable containing the domain-specific sampling variances of the direct estimates that are included in dt
#' @param transformation character, either "no" (default) or "arcsin".
#' @param eff_smpsize character, name of the variable containing the effective sample sizes that are included in dt. Required argument when the arcsin transformation is chosen. Defaults to NULL.
#' 
#' @import data.table
#' @import caret
#' @importFrom emdi step fh

step_wrapper_fh <- function(dt, xvars, y, cor_thresh = 0.95, criteria = "AIC",
                         vardir, transformation = "no", eff_smpsize) {
  
  dt <- as.data.table(dt)
  
  # Drop columns that are entirely NA
  dt <- dt[, which(unlist(lapply(dt, function(x) !all(is.na(x))))), with = FALSE]
  
  xvars <- xvars[xvars %in% colnames(dt)]
  
  # Keep only complete cases
  dt <- dt[complete.cases(dt),] 
  
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
  
  # Stepwise selection
  fh_args <- list(
    fixed = model_formula,
    vardir = vardir,
    combined_data = dt,
    method = "ml",
    MSE = FALSE
  )
  
  if (transformation == "arcsin") {
    fh_args$transformation <- "arcsin"
    fh_args$backtransformation <- "bc"
    fh_args$eff_smpsize <- eff_smpsize
  } else {
    fh_args$transformation <- "no"
    fh_args$B <- c(0, 50)
  }
  
  fh_model <- do.call(emdi::fh, fh_args)
  
  stepwise_model <- emdi::step(fh_model, criteria = criteria)
  
  return(stepwise_model)
  
}
