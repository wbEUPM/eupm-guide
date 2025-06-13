outdetect <- function(x, verbose = TRUE) {
  if (!requireNamespace("car", quietly = TRUE)) stop("Package 'car' is required.")
  if (!requireNamespace("robustbase", quietly = TRUE)) stop("Package 'robustbase' is required.")
  if (!requireNamespace("bestNormalize", quietly = TRUE)) stop("Package 'bestNormalize' is required.")

  # Check input
  if (all(is.na(x))) stop("Input contains only NA values.")
  
  nonzero_vals <- x[x != 0 & !is.na(x)]
  if (length(nonzero_vals) < 1) stop("No non-zero observations available.")
  
  # Estimate λ and transform non-zero values
  pt <- car::powerTransform(nonzero_vals ~ 1, family = "yjPower")
  lambda <- pt$lambda
  x_nonzero_YJ <- car::yjPower(nonzero_vals, lambda = lambda)
  
  # Compute robust location and scale
  med_YJ <- median(x_nonzero_YJ)
  qn <- robustbase::Qn(x_nonzero_YJ, constant = 2.2219, finite.corr = FALSE)
  if (verbose) message(sprintf("Qn value: %.4f", qn))
  
  # Apply same transformation to all x
  x_all_YJ <- car::yjPower(x, lambda = lambda)
  
  # Compute z-scores and flag outliers
  z_scores <- (x_all_YJ - med_YJ) / qn
  outlier_flag <- rep(0, length(x))
  outlier_flag[!is.na(z_scores) & round(z_scores, 2) >= 2.99] <- 2
  outlier_flag[!is.na(z_scores) & round(z_scores, 2) <= -2.99] <- 1
  outlier_flag[x == 0 & !is.na(x)] <- 0
  
  # Return data frame
  return(data.frame(
    value = x,
    z = ifelse(is.na(x), NA, z_scores),
    outlier = outlier_flag
  ))
}

