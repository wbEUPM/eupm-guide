# ------------------------------------------------------------------------------
# Function: mg
# Purpose : Merge two datasets, keeping master version of overlapping variables
#
# Revision History:
# Date       | Author   | Description
#------------|----------|-----------------------------------------------
# 2024-03-12 | ST       | Initial version with summary table and merge codes
# 2025-03-20 | Rose     | ?
# 2025-03-26 | ST       | 1.Modified to retain master_data columns on conflict as 
#            |          | default.   
#            |          | 2. add the option to keep both using x.var & y.var1.   
#            |          | 3. Replace var merge when var merge already exists.   
# 2025-05-16 | ST       | Added lines to install missing packages. 
# 2025-06-09 | ST       | Improved the lines to install missing packages. 
# 2025-06-11 | ST       | Added to option to drop merge 
# ------------------------------------------------------------------------------

ensure_packages <- function(packages) {
  to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(to_install) > 0) install.packages(to_install)
  invisible(lapply(packages, library, character.only = TRUE))
}

ensure_packages(c("dplyr", "rlang", "knitr"))




mg <- function(master_data, using_data, key, keep_both = FALSE, drop_merge = FALSE) {
  
  # Handle potential conflict with existing "merge" column
  if ("merge" %in% names(master_data)) {
    warning("'merge' column already exists in master_data and will be replaced.")
    master_data <- master_data %>% select(-merge)
  }
  if ("merge" %in% names(using_data)) {
    warning("'merge' column already exists in using_data and will be replaced.")
    using_data <- using_data %>% select(-merge)
  }
  
  # Add identifiers for tracking match status
  master_data <- master_data %>% mutate(master_v = 1)
  using_data  <- using_data %>% mutate(using_v = 1)
  
  # Identify common columns (excluding keys)
  common_cols <- setdiff(intersect(names(master_data), names(using_data)), key)
  
  # Drop common columns from using_data if keep_both is FALSE
  if (!keep_both && length(common_cols) > 0) {
    using_data <- using_data %>% select(-all_of(common_cols))
  }
  
  # Perform full join and generate new "merge" status column
  merged_data <- master_data %>%
    full_join(using_data, by = key) %>%
    mutate(merge = case_when(
      is.na(using_v)  ~ 1,
      is.na(master_v) ~ 2,
      TRUE            ~ 3
    )) %>%
    select(-master_v, -using_v)
  
  # Show merge summary
  merge_summary <- merged_data %>%
    group_by(merge) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    mutate(Result = case_when(
      merge == 1 ~ "Not matched from master",
      merge == 2 ~ "Not matched from using",
      merge == 3 ~ "Matched",
      TRUE       ~ "Unknown"
    ))
  
  cat("Merge result:\n")
  print(kable(merge_summary))
  
  # Drop merge column if requested
  if (drop_merge) {
    merged_data <- merged_data %>% select(-merge)
  }
  
  return(merged_data)
}


#Example Data
# After saving this file to your work directory, include the following line:
# df1 <- data.frame(id = c(1, 2, 3, 4),
#                   name = c("Alice", "Bob", "Charlie", "David"),
#                   dup = c(1,1,1,1))
# df2 <- data.frame(id = c(1, 2, 3),
#                   age = c(25, 30, 35),
#                   dup = c(2,2,2),
#                   merge = c(1,1,1))
# 
# 
# # Ver 1: Run the Function (Default)
# merged <- mg(df1, df2, key="id")
# # View Merged Data
# print(merged)
# # Confirm that dup from using_data is not included.

# Drop unmatched observations from the master dataset (merge != 1)
# merged_cleaned <- merged %>% filter(merge != 1)
# 
# # View Cleaned Data
# print(merged_cleaned)
# 
# # Drop merge result variable
# merged_cleaned2 <- merged %>% select(-merge)
# # View Cleaned Data 2
# print(merged_cleaned2)
# 
# 
# 
# # Ver 2: Run the Function using keep_both option
# merged <- mg(df1, df2, key="id", keep_both = TRUE)
# # View Merged Data
# print(merged)
# # Confirm that dup from using_data is not included.
# 
# # Drop unmatched observations from the master dataset (merge != 1)
# merged_cleaned <- merged %>% filter(merge != 1)
# 
# # View Cleaned Data
# print(merged_cleaned)
# 
# # Drop merge result variable
# merged_cleaned2 <- merged %>% select(-merge)
# # View Cleaned Data 2
# print(merged_cleaned2)