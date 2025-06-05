################################################################################
############## PREPARE THE NGA GHS DATA TO SHOW MFH MODELLING ##################
################################################################################

pacman::p_load(tidyverse, sf, data.table, haven)

### ok lets start with the 2018 data

ghsp4_dt <- haven::read_dta("data/nga/NGA_2018_GHSP-W4_v01_M_v01_A_SSAPOV_GMD.dta")

ghsp3_list <- 
  list.files(path = "data/nga",
             pattern = "^NGA_2015_GHSP-W3_v01_M_v01_.*\\.dta$",
             full.names = TRUE) %>% 
  lapply(X = .,
         FUN = haven::read_dta)

totcons4_dt <- haven::read_dta("data/nga/totcons_final.dta")

labor4_dt <- haven::read_dta("data/nga/sect3a_harvestw4.dta")
  

### prepare names for the ghsp3 data
module_list <- 
  list.files(path = "data/nga",
             pattern = "^NGA_2015_GHSP-W3_v01_M_v01_.*\\.dta$") %>%
  gsub(pattern = "NGA_2015_GHSP-W3_v01_M_v01_A_",
       replacement = "",
       x = .) %>%
  gsub(pattern = ".dta",
       replacement = "",
       x = .)

names(ghsp3_list) <- module_list

### lets combine all the data with a function

combine_gmd <- function(dt_list){
  
  dt <- 
    dt_list$SSAPOV_I %>%
    merge(dt_list$SSAPOV_L, all = TRUE) %>%
    merge(dt_list$SSAPOV_P, all = TRUE) %>%
    merge(dt_list$SSAPOV_H, all = TRUE)

  return(dt)    
  
}

### put together the wave 3 data
ghsp3_dt <- combine_gmd(ghsp3_list)


### lets do the same for wave 2
ghsp2_list <- 
  list.files(path = "data/nga",
             pattern = "^NGA_2012_GHSP-W2_v01_M_v01_.*\\.dta$",
             full.names = TRUE) %>% 
  lapply(X = .,
         FUN = haven::read_dta)

module_list <- 
  list.files(path = "data/nga",
             pattern = "^NGA_2012_GHSP-W2_v01_M_v01_.*\\.dta$") %>%
  gsub(pattern = "NGA_2012_GHSP-W2_v01_M_v01_A_",
       replacement = "",
       x = .) %>%
  gsub(pattern = ".dta",
       replacement = "",
       x = .)

names(ghsp2_list) <- module_list

ghsp2_dt <- combine_gmd(ghsp2_list)


### include poverty and geocode data thats missing from GHSP4

ghsp4_dt <- merge(ghsp4_dt, totcons4_dt)


ghsp4_dt <- 
  ghsp4_dt %>%
  rename(wel_PPPnom = "totcons_pc") %>%
  rename(region2 = "subnatid2") %>%
  rename(region1 = "subnatid1") %>%
  rename(cluster = "ea") |>
  rename(relathh6 = relationharm)

ghsp2_dt <- 
  ghsp2_dt |>
  mutate(age = ageyrs)

ghsp3_dt <- 
  ghsp3_dt |>
  mutate(age = ageyrs)

### include cpi and ppp values for each survey
ghsp4_dt <- 
  ghsp4_dt %>%
  mutate(cpi = 1.209061823059,
         ppp = 123.67238)

ghsp3_dt <- 
  ghsp3_dt %>%
  mutate(cpi = 1.209061823059,
         ppp = 123.67238)

ghsp2_dt <- 
  ghsp2_dt %>%
  mutate(cpi = 1.209061823059,
         ppp = 123.67238)


ghsp4_dt <- 
  ghsp4_dt %>%
  rename(wta_hh = "wt_wave4")

comnames_list <- Reduce(f = intersect,
                        x = list(colnames(ghsp2_dt),
                                 colnames(ghsp3_dt),
                                 colnames(ghsp4_dt)))


#### ok lets combine the data and create some variables

find_inconsistent_datatypes <- function(dt_list, varnames) {
  # Get class info for each dataset
  class_list <- lapply(dt_list, function(df) {
    sapply(df[, varnames, drop = FALSE], function(x) paste(class(x), collapse = "/"))
  })
  
  # Convert to data.frame with variables as rows and datasets as columns
  class_df <- as.data.frame(do.call(cbind, class_list), stringsAsFactors = FALSE)
  
  # Add variable names as a column
  class_df$varname <- varnames
  
  # Reorder columns for readability
  class_df <- class_df[, c(ncol(class_df), 1:(ncol(class_df) - 1))]
  
  # Name columns with dataset names, if provided
  if (!is.null(names(dt_list))) {
    colnames(class_df)[-1] <- names(dt_list)
  } else {
    colnames(class_df)[-1] <- paste0("df", seq_along(dt_list))
  }
  
  # Identify inconsistent types
  is_inconsistent <- apply(class_df[, -1], 1, function(x) length(unique(x)) > 1)
  inconsistent_df <- class_df[is_inconsistent, ]
  
  return(inconsistent_df)
}

drop_vars <- 
find_inconsistent_datatypes(dt_list = list(ghsp2 = ghsp2_dt, 
                                           ghsp3 = ghsp3_dt, 
                                           ghsp4 = ghsp4_dt),
                            varnames = comnames_list) |>
  dplyr::select(varname) |>
  as.list() |>
  unlist() |>
  unname()

comnames_list <- comnames_list[!comnames_list %in% drop_vars]


ghsp2_dt_clean <- ghsp2_dt |> mutate(across(all_of(comnames_list), zap_labels))
ghsp3_dt_clean <- ghsp3_dt |> mutate(across(all_of(comnames_list), zap_labels))
ghsp4_dt_clean <- ghsp4_dt |> mutate(across(all_of(comnames_list), zap_labels))

ghsp_dt <- bind_rows(
  ghsp2_dt_clean |> select(all_of(comnames_list)),
  ghsp3_dt_clean |> select(all_of(comnames_list)),
  ghsp4_dt_clean |> select(all_of(comnames_list))
)

rm(ghsp2_dt_clean, ghsp3_dt_clean, ghsp4_dt_clean)



#### lets prepare variables for small area estimation

### start with the age category variables

ghsp_dt <- 
  ghsp_dt %>%
  mutate(agecat = gsub(pattern = "-", "", agecat),
         agecat = gsub(pattern = " ", "", agecat))
  
ghsp_dt <- 
  ghsp_dt %>%
  mutate(agecat = ifelse(agecat == "" | is.na(agecat), "Notstated", agecat))


### the variables we have (age)













