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
  vcov_list <- lapply(X = dom_list, FUN = surv_vcov)
  var_dt <- Reduce(f = "rbind", x = vcov_list) |> 
    as_tibble() |> 
    mutate({{area_var}} := dom_list)
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
