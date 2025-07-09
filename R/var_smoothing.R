#' A function to perform variance smoothing
#' 
#' The variance smoothing function applies the methodology of You and Hiridoglou (2023)
#' which uses simply log linear regression to estimate direct variances for sample 
#' poverty rates which is useful for replacing poverty rates in areas with low sampling.
#' 
#' @param domain a vector of unique domain/target areas
#' @param direct_var the raw variances estimated from sample data
#' @param sampsize the sample size for each domain
#' 
#' @export

varsmoothie_king <- function(domain,
                             direct_var,
                             sampsize){
  
  dt <- data.table(Domain = domain,
                   var = direct_var,
                   n = sampsize)
  
  dt$log_n <- log(dt$n)
  dt$log_var <- log(dt$var)
  
  lm_model <- lm(formula = log_var ~ log_n,
                 data = dt[!(abs(dt$log_var) == Inf),])
  
  dt$pred_var <- predict(lm_model, newdata = dt)
  residual_var <-  summary(lm_model)$sigma^2
  dt$var_smooth <- exp(dt$pred_var) * exp(residual_var/2)
  
  return(dt[, c("Domain", "var_smooth"), with = F])
  
}