################################################################################
########### ESTIMATING THE FH MODEL WITH SPAIN DATA ###########
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

# Aggregate auxiliary variables
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
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)),
            N = n()) %>%
  ungroup()
#weighted.mean(.x, na.rm = TRUE, w = weight)),


############### Poverty rate

#### step 1: compute the direct estimates for the poverty rates

income_dt$poor <- as.numeric(as.integer(income_dt$income2012 < income_dt$povline2012))
table(income_dt$poor)

svy_design <- svydesign(
  id = ~1,
  weights = ~weight,
  #fpc = ~N,
  data = income_dt
)

overallmean <- svymean(~income_dt$poor, design = svy_design, deff = TRUE)
# mean 0.2165049


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
  left_join(as.data.frame(table(income_dt$prov)) %>%
              rename(prov = Var1) %>%
              mutate(prov = as.integer(prov)) %>%
              rename(n = Freq), by = "prov")


comb_Data_poor <- dir_poor %>%
  left_join(aux_agg, by = "prov")

#comb_Data <- combine_data(pop_data = aux_agg, pop_domains = "prov", 
 #                         smp_data = dir_poor, smp_domains = "prov")

saveRDS(comb_Data_poor, file = "data/comb_Data_poor.RDS")

################################################################################
### Estimation of FH model

comb_Data_poor <- readRDS("data/comb_Data_poor.RDS")

comb_Data_poor <- comb_Data_poor %>%
  mutate(n_eff = n/DEff.poor)

fh_start <- step(fh(
  fixed = poor ~ 
    #gen + age2 + 
  #  age3 + age4 + age5 +
   # + educ1 + 
   # educ2 + educ3 +
  #  nat1 + 
  #  labor1 +
  #  labor2,
    abs + ntl + aec + schyrs + mkt,
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
estimators(fh_arcsin, MSE = TRUE, CV = TRUE)

data("sizeprov")
comb_Data_poor$ratio_n <- sizeprov$Nd/(sum(sizeprov$Nd))


fh_bench <- benchmark(fh_arcsin,
                      benchmark = 0.2165049,
                      share = comb_Data_poor$ratio_n, 
                      type = "ratio",
                      overwrite = TRUE)

saveRDS(fh_arcsin, "data/fh_arcsin.RDS")

############### Mean income

emdi_direct <- direct(
  y = "income2012", smp_data = income_dt, smp_domains = "prov", 
  weights = "weight", #threshold = 11064.82,
  var = TRUE, boot_type = "naive", B = 50, seed = 123, X_calib = NULL,
  totals = NULL, na.rm = TRUE
)

dir_mean <- emdi_direct$ind %>%
  dplyr::select(Domain, Mean) %>%
  left_join(emdi_direct$MSE %>%
              dplyr::select(Domain, Mean) %>%
              rename(vardir = Mean), by = "Domain") %>%
  rename(prov = Domain) %>%
  mutate(prov = as.integer(prov))


domsizeMean <- income_dt %>% group_by(prov) %>%
  summarise(sumW = sum(weight)) %>%
  ungroup() %>% data.frame()

surveyMean <- sae::direct(y = as.numeric(income_dt$income2012),
                          dom = factor(income_dt$prov),
                          sweight = income_dt$weight,
                          domsize = domsizeMean)

surveyMean <- surveyMean %>%
  rename(prov = Domain, n = SampSize, Mean = Direct,
         mean_sd = SD, mean_cv = CV) %>%
  mutate(mean_cv = mean_cv / 100) %>%
  mutate(vardir = mean_sd^2) %>%
  mutate(prov = as.integer(prov))

comb_Data_mean <- surveyMean %>%
  left_join(aux_agg, by = "prov")

#comb_Data <- combine_data(pop_data = aux_agg, pop_domains = "prov", 
#                         smp_data = dir_poor, smp_domains = "prov")

saveRDS(comb_Data_mean, file = "data/comb_Data_mean.RDS")

################################################################################
### Estimation of FH model

comb_Data_mean <- readRDS("data/comb_Data_mean.RDS")

fh_start <- step(fh(
  fixed = Mean ~ 
    gen + age2 + 
    age3 + age4 + age5 +
    + educ1 + 
    educ2 + educ3 +
    nat1 + 
    labor1 +
    labor2 +
    abs + ntl + aec + schyrs + mkt,
  vardir = "vardir", combined_data = comb_Data_mean, domains = "prov",
  method = "ml", transformation = "log", backtransformation = "bc_sm",
  MSE = FALSE)) 

fh_log <- fh(
  fixed = formula(fh_start),
  vardir = "vardir", combined_data = comb_Data_mean, domains = "prov",
  method = "ml", transformation = "log", backtransformation = "bc_sm",
  MSE = TRUE)

fh_log$model$variance
summary(fh_log$ind)
summary(fh_log$MSE)
summary(fh_log)
plot(fh_log)
compare(fh_log)
compare_plot(fh_log, MSE = TRUE, CV = TRUE)
estimators(fh_log, MSE = TRUE, CV = TRUE)

data("sizeprov")
comb_Data_mean$ratio_n <- sizeprov$Nd/(sum(sizeprov$Nd))

income_dt <- readRDS("data/incomedata.RDS")
mean(income_dt$income2012)
# 12109.51

fh_bench <- benchmark(fh_log,
                      benchmark = 12109.51,
                      share = comb_Data_mean$ratio_n, 
                      type = "ratio",
                      overwrite = TRUE)

saveRDS(fh_log, "data/fh_log.RDS")
