################################################################################
########### ESTIMATING THE MFH SPATIO TEMPORAL MODEL WITH SPAIN DATA ###########
################################################################################

### load the libraries

if (sum(installed.packages()[,1] %in% "pacman") != 1){
  
  install.packages("pacman")
  
}

pacman::p_load(sf, data.table, tidyverse, car, msae, sae, survey, spdep, emdi)


### read in the datasets

income_dt <- readRDS("data/incomedata.RDS")

#shp_dt <- readRDS("data/shapes/spainshape.RDS")


##### let us start by developing the typical FH model 

#### step 1: compute the direct estimates for the poverty rates

# Draw 25% sample per province from incomedata
n_per_province <- round(0.25*(as.data.frame(table(income_dt$prov))$Freq), 0) # adjust as needed
set.seed(123)
sample_id <- stratsrs(income_dt$prov, n_per_province) 
samp_data <- income_dt[sample_id,]

samp_data$poor <- as.numeric(as.integer(samp_data$income2012 < samp_data$povline2012))
table(samp_data$poor)

svy_design <- svydesign(
  id = ~1,
  weights = ~weight,
  #fpc = ~N,
  data = samp_data
)

overallmean <- svymean(~samp_data$poor, design = svy_design, deff = TRUE)

result <- svyby(
  ~poor,
  by = ~prov,
  design = svy_design,
  FUN = svymean,
  na.rm = TRUE, 
  deff = TRUE
)
result
dir_poor <- result %>%
  mutate(vardir = se^2) %>%
  left_join(as.data.frame(table(samp_data$prov)) %>%
              rename(prov = Var1) %>%
              mutate(prov = as.integer(prov)) %>%
              rename(n = Freq), by = "prov")

colnames(income_dt)
summary(income_dt)

### let create some additional variables
# change dummies of gen and nat to 0 and 1
income_dt <- 
  income_dt %>%
  mutate(across(c(gen, nat), ~ case_when(
    .x == 1 ~ 0,
    .x == 2 ~ 1,
    TRUE ~ NA_real_
  )))


### create the province level dataset

candidate_vars <- colnames(income_dt)[!colnames(income_dt) %in% 
                                        c("provlab", "prov", "income",
                                          "weight", "povline", "y",
                                          "poverty", "y0", "y1", "y2",
                                          "p0_prov", "p1_prov", "p2_prov",
                                          "ac", "nat", "educ", "labor",
                                          "age", "poor")]

candidate_vars <- candidate_vars[!grepl("sampsize|prov|poverty|income|povline|^v[0-9]|^y[0-9]", 
                                        candidate_vars)]

aux_agg <- income_dt %>%
  dplyr::select(weight, prov, all_of(candidate_vars)) %>%
  group_by(prov) %>%
  summarise(across(everything(), ~ weighted.mean(.x, na.rm = TRUE, w = weight)),
            N = n()) %>%
  ungroup()

comb_Data_poor <- dir_poor %>%
  left_join(aux_agg, by = "prov")

#comb_Data <- combine_data(pop_data = aux_agg, pop_domains = "prov", 
 #                         smp_data = dir_poor, smp_domains = "prov")

saveRDS(comb_Data_poor, file = "data/shapes/comb_Data_poor.RDS")

################################################################################
### Estimation of FH model

comb_Data_poor <- readRDS("data/shapes/comb_Data_poor.RDS")

comb_Data_poor <- comb_Data_poor %>%
  mutate(n_eff = n/DEff.poor)

fh_start <- step(fh(
  fixed = poor ~ gen + #age2 + 
    age3 + age4 + age5 +
   # + educ1 + 
    educ2 + educ3 +
    nat1 + 
    #labor1 +
    labor2,
    #abs + ntl + aec + schyrs + mkt,
  vardir = "vardir", combined_data = comb_Data_poor, domains = "prov",
  method = "ml", transformation = "arcsin", backtransformation = "bc",
  eff_smpsize = "n_eff", MSE = FALSE)) 

fh_arcsin <- fh(
  fixed = formula(fh_start),
  vardir = "vardir", combined_data = comb_Data_poor, domains = "prov",
  method = "ml", transformation = "arcsin", backtransformation = "bc",
  eff_smpsize = "n_eff", MSE = TRUE, mse_type = "boot", B = c(50, 0))

summary(fh_arcsin$ind)
summary(fh_arcsin$MSE)
summary(fh_arcsin)
plot(fh_arcsin)
compare(fh_arcsin)
compare_plot(fh_arcsin, MSE = TRUE, CV = TRUE)
saveRDS(model3_obj, "data/modelmfh3.RDS")
