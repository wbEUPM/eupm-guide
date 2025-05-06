################################################################################
################# PREPARE THE ALBANIA DATA FOR MFH MODELLING ###################
################################################################################

pacman::p_load(haven, data.table, tidyverse)

alb_dt <- haven::read_dta("data/alb_2005_2008.dta")

### quickly write the labels to file 




