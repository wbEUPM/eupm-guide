#' Convert text into expression (helper)
#' @description
#' A little worker function to convert expressions that are character to the 
#' right format for the survey package
#' 
#' @importFrom rlang expr
convert_expr <- function(x){
  y <- if (is.character(x)){
    rlang::expr(~!!sym(x))
  } else {
    rlang::expr(~ !!x)
  }
  y <- eval(y)
  return(y)
}


#' Estimate variance covariance matrix for Multivariate Fay Herriot Model
#' 
#' A function to estimate the variance covariance matrix for the multivariate fay herriot
#' model 
#' 
#' @param dt, `data.frame` the individual/household unit level data with the variable of interest from
#' which the covariance matrix will be computed
#' @param domain `character` column name for the domain variable in `dt`
#' @param ids `character` the cluster variable
#' @param weights `character` the weight variable
#' @param fpc the finite population correction variable, see `survey::svymean()` for more details
#' @param strata `character` strata variable
#' @param probs see `survey::svymean()` for more details
#' @param yvars the outcome variable of interest
#' @param deff `logical` for whether or not to estimate the design effect
#' 
#' @export 
#' 
#' @import survey
#' 

compute_vcov <- function(dt,
                         domain,
                         ids,
                         weights = NULL,
                         fpc = NULL,
                         strata = NULL,
                         probs = NULL,
                         yvars,
                         deff = FALSE,
                         ...){
  
  
  ### run it on the proper arguments
  ids <- convert_expr(x = ids)
  if (is.null(weights)){
    dt$weights <- 1
  }
  weights <- convert_expr(x = weights)
  if (!is.null(fpc)){    
    fpc <- convert_expr(x = fpc)
  }
  if (!is.null(probs)){
    probs <- convert_expr(x = probs)
  }
  if (!is.null(strata)){
    strata <- convert_expr(x = strata)
  }
  
  ### create a survey design object from the survey package 
  dom_list <- unique(dt[[domain]])
  
  surv_vcov <- function(x){
    
    design_obj <- svydesign(ids = ids,
                            probs = probs,
                            strata = strata,
                            fpc = fpc,
                            data = dt |>
                              dplyr::filter(!!sym(domain) == x),
                            weights = weights)
    
    ### prepare the y formula for the svymean function
    yform <- reformulate(yvars)
    var_dt <- svymean(yform, design_obj, ...) ## use the svymean object to compute variance
    var_dt <- vcov(var_dt) ### compute variance covariance matrix
    var_dt <- as.numeric(c(diag(var_dt), var_dt[lower.tri(var_dt, diag = FALSE)]))
    pair_list <- c(lapply(yvars, rep, 2), combn(yvars, 2, simplify = FALSE))
    pair_list <- lapply(pair_list,
                        function(x){
                          zz <- paste0("v_", paste(x, collapse = "__"))
                          return(zz)
                        })
    pair_list <- unlist(pair_list)
    names(var_dt) <- pair_list
    return(var_dt)
  }
  var_dt <- 
    lapply(X = dom_list, FUN = surv_vcov) |>
    lapply(FUN = function(x){
      
      y <- x |> t() |> as_tibble()
      
      return(y)
      
    }) |>
    Reduce(f = rbind) |>
    mutate(!!sym(domain) := dom_list) |>
    dplyr::select(!!sym(domain), starts_with("v_"))
  
  return(var_dt)  
}


#' Estimate variance covariance matrix for pairs of observations
compute_vcov_pairs <- function(dt,
                         domain,
                         ids,
                         weights = NULL,
                         fpc = NULL,
                         strata = NULL,
                         probs = NULL,
                         yvars,
                         deff = FALSE,
                         ...){
  
  
  ### run it on the proper arguments
  ids <- convert_expr(x = ids)
  if (is.null(weights)){
    dt$weights <- 1
  }
  weights <- convert_expr(x = weights)
  if (!is.null(fpc)){    
    fpc <- convert_expr(x = fpc)
  }
  if (!is.null(probs)){
    probs <- convert_expr(x = probs)
  }
  if (!is.null(strata)){
    strata <- convert_expr(x = strata)
  }
  
  ### create a survey design object from the survey package 
  dom_list <- unique(dt[[domain]])
  all_pairs <- combn(yvars, 2, simplify = FALSE)
  
  surv_vcov_2 <- function(dta, yvars_local, ...){
    yform <- reformulate(yvars_local)
    var_covar_name <- paste0("v_", paste(yvars_local, collapse = "__"))
    var_n_name <- paste0("n_", paste(yvars_local, collapse = "__"))
    local_dta <- dta |> filter(if_all(all_of(yvars_local), ~!is.na(.)))
    local_n_obs <- nrow(local_dta)
    out <- tibble(x = 0, n = local_n_obs) 
    if (local_n_obs > 1) {
      design_obj <- svydesign(
        ids = ids,
        variables = yform,
        probs = probs,
        strata = strata,
        fpc = fpc,
        data = local_dta,
        weights = weights
      )
      var_mean <- svymean(yform, design_obj, ...) 
      # var_var <- svyvar(yform, design_obj, ...) 
      var_covar <- vcov(var_mean)
      out <- 
        tibble(x = var_covar[yvars_local[[1]], yvars_local[[2]]], n = local_n_obs)
    }
    out |> rename({{var_covar_name}} := x, {{var_n_name}} := n)
  }
  dt |>
    group_by(pick(all_of(domain))) |> 
    nest() |> 
    mutate(
      data = map(data, ~{
        # browser()
        dta_local <- .x
        all_pairs |>
          map(~{surv_vcov_2(dta = dta_local, yvars_local = .x, na.rm = TRUE)}) |> 
          bind_cols()
      })
    ) |> 
    bind_rows() |> 
    unnest(data)
}



#' Generate Domain-Level Variance-Covariance Matrix
#'
#' Computes a variance-covariance matrix of direct survey estimates at the domain level, with optional 
#' functionality to smooth variances, replace outliers, or impose a constant covariance structure.
#'
#' @param dt A data frame or tibble containing the microdata used to compute direct estimates.
#' @param domain A character string specifying the column name for domain/group identifiers.
#' @param ids Optional. Variable name for primary sampling unit (PSU) identifiers. Default is `NULL`.
#' @param weights Optional. Variable name for survey weights. Default is `NULL`.
#' @param fpc Optional. Variable name for finite population correction (FPC). Default is `NULL`.
#' @param strata Optional. Variable name for survey strata. Default is `NULL`.
#' @param probs Optional. Not currently used.
#' @param yvars A character vector of outcome variable names to include in the variance-covariance matrix.
#' @param deff Logical. If `TRUE`, adjusts variances using design effect estimation. Default is `FALSE`.
#' @param constant_cov Logical. If `TRUE`, replaces off-diagonal covariances with a constant value. Default is `FALSE`.
#' @param cov_value Numeric. Value to assign to off-diagonal covariances if `constant_cov = TRUE`.
#' @param use_smooth_var Logical. If `TRUE`, returns only smoothed variances using the `varsmoothie_king()` function. 
#' Cannot be `TRUE` if `replace_outlier = TRUE`. Default is `FALSE`.
#' @param replace_outlier Logical. If `TRUE`, replaces variance outliers (based on IQR rule) 
#' with their smoothed counterparts. Cannot be `TRUE` if `use_smooth_var = TRUE`. Default is `FALSE`.
#'
#' @return A tibble containing either:
#' - Domain-level variance-covariance matrix (`v_` columns), optionally with outliers replaced by smoothed estimates
#' - Or a smoothed matrix (`vsv_` columns) if `use_smooth_var = TRUE`
#'
#' @details
#' The function uses the `compute_vcov()` helper to compute design-based variance and covariance estimates 
#' for all combinations of `yvars`. Sample size per domain is used for smoothing. If `replace_outlier = TRUE`, 
#' outliers beyond the standard IQR rule are substituted with their smoothed versions. If `use_smooth_var = TRUE`, 
#' the function returns only smoothed variances.
#'
#' @seealso [`compute_vcov()`], [`varsmoothie_king()`]
#'
#' @export
#' 
#' @examples
#' 
#' ## compute the standard variance covariance matrix
#' 
#' time <- 2012:2014
#' 
#' 
#' income_dt <- reduce(
#'   time,
#'   .init = income_dt,
#'   .f = function(df, y) {
#'     income_var <- sym(paste0("income", y))
#'     povline_var <- sym(paste0("povline", y))
#'     new_var <- paste0("poor", y)
#'     
#'     df %>%
#'       mutate(!!new_var := if_else(!!income_var < !!povline_var, 1, 0))
#'   }
#' )
#' 
#' ## compute standard variance-covariance matrix
#' gen_vcov(dt = income_dt, domain = "prov", ids = "psu", weights = "weight", 
#'          yvars = c("poor2012", "poor2013", "poor2014"))
#' 
#' ## replace covriance portion with constant fixed value
#' gen_vcov(dt = income_dt, domain = "prov", ids = "psu", weights = "weight", 
#'          yvars = c("poor2012", "poor2013", "poor2014"),
#'          constant_cov = TRUE, cov_value = 0)
#'          
#' ## replace outliers in variance-covariance matrix with smoothed values
#' gen_vcov(dt = income_dt, domain = "prov", ids = "psu", weights = "weight", 
#'          yvars = c("poor2012", "poor2013", "poor2014"),
#'          replace_outlier = TRUE)
#'          
#' ## replace ALL variance-covariance matrix entirely with smoothed values
#' gen_vcov(dt = income_dt, domain = "prov", ids = "psu", weights = "weight", 
#'          yvars = c("poor2012", "poor2013", "poor2014"),
#'          use_smooth_var = TRUE)
#' 
#' 
#' 
#' 
#' 




gen_vcov <- function(dt,
                     domain,
                     ids = NULL,
                     weights = NULL,
                     fpc = NULL,
                     strata = NULL,
                     probs = NULL,
                     yvars,
                     deff = FALSE,
                     constant_cov = FALSE,
                     cov_value,
                     use_smooth_var = FALSE,
                     replace_outlier = FALSE){
  
  if (replace_outlier && use_smooth_var) {
    stop("You cannot set both `replace_outlier = TRUE` and `use_smooth_var = TRUE`. Please choose only one.")
  }
  
  ### compute variance covariance matrix
  vcov_dt <- compute_vcov(dt = dt,
                          domain = domain,
                          ids = ids,
                          weights = weights,
                          fpc = fpc,
                          strata = strata,
                          probs = probs,
                          yvars = yvars,
                          deff = deff)
  
  ### merge in sample size
  vcov_dt <- left_join(vcov_dt, 
                       dt |> 
                         count(!!sym(domain), 
                               name = "sampsize"), 
                       by = domain)
  
  if (constant_cov == TRUE){
    
    vars <- 
      colnames(vcov_dt)[grepl("^v_", colnames(vcov_dt))] %>%
      tibble(name = gsub("^v_", "", x = .)) %>%
      separate(name, into = c("var1", "var2"), sep = "__") |>
      rename(var = ".") |>
      mutate(equal = ifelse(var1 == var2, 1, 0)) |>
      filter(equal == 0) |>
      dplyr::select(var) |>
      c() |>
      unlist()
    
    vcov_dt[, vars] <- cov_value
    
  }
  
  
  smooth_vcov <- function(){
    
    vardir <- grep("^v_", names(vcov_dt), value = TRUE)
    
    vcov_dt <-
      lapply(X = vardir,
             FUN = function(x){
               
               z <- varsmoothie_king(domain = vcov_dt[[domain]],
                                     direct_var = vcov_dt[[x]],
                                     sampsize = vcov_dt[["sampsize"]]) |>
                 as.data.table() |>
                 setnames(old = "var_smooth", new = paste0("vs", x)) |>
                 as_tibble()
               
               
               return(z)
               
             }) %>%
      Reduce(f = "merge",
             x = .) %>%
      merge(x = vcov_dt,
            y = .,
            by.y = "Domain",
            by.x = domain) |>
      as_tibble()
    
    return(vcov_dt)
    
  }
  
  
  
  if (replace_outlier == TRUE){
    
    message("Replacing outliers with smoothed estimates...")
    
    vcov_dt <- smooth_vcov()
    
    vardir <- grep("^v_", names(vcov_dt), value = TRUE)
    
    iqr_bounds <- 
      map(vcov_dt[vardir], 
          ~ {
            q1 <- quantile(.x, 0.25, na.rm = TRUE)
            q3 <- quantile(.x, 0.75, na.rm = TRUE)
            iqr <- q3 - q1
            list(lower = q1 - 1.5 * iqr, upper = q3 + 1.5 * iqr)
          }
      )
    
    # Step 3: Create named list for quick lookup
    names(iqr_bounds) <- vardir
    
    # Step 4: Replace outliers using mutate + across
    vcov_dt <- 
      vcov_dt %>%
      mutate(across(
        all_of(vardir),
        ~ {
          bounds <- iqr_bounds[[cur_column()]]
          smoothed_col <- get(paste0("vsv_", str_remove(cur_column(), "^v_")))
          if_else(.x < bounds$lower | .x > bounds$upper, smoothed_col, .x)
        },
        .names = "{.col}"
      )) |>
      dplyr::select(!!sym(domain), starts_with("v_"))
    
    
    
  }
  
  
  if (use_smooth_var == TRUE){
    
    message("Using fully smoothed variance-covariance matrix...")
    
    vcov_dt <- smooth_vcov()
    
    vcov_dt <- 
      vcov_dt |>
      dplyr::select(!!sym(domain), starts_with("vsv_"))
    
  }
  
  return(vcov_dt)
  
}






















