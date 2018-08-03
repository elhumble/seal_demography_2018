# collect best lhoods from different fastsimcoal runs

library(readr)
library(purrr)
library(stringr)
library(parallel)
options(scipen=999)
setwd(paste0("/home/emily/angsd_pilon/"))
# how many bootstraps
# nboot <- 100
# num_sim <- 50


#~~ Bottleneck model

# check folders in fsc_run/
folders <- list.files("fsc_run/", pattern = "fscrun[0-9]")

# collect likelihoods 
collect_vars <- function(folder){
  fsc_pars <- read_delim(paste0("fsc_run/", folder, "/afs/afs.bestlhoods"), delim = "\t")
}

all_vars <- folders %>% 
  map(safely(collect_vars)) %>% 
  purrr::transpose()

all_vars_working <- do.call(rbind, all_vars[[1]]) %>% 
  dplyr::mutate(lh_diff = MaxEstLhood - MaxObsLhood)

all_vars_working[which.max(all_vars_working$MaxEstLhood), ]

hist(all_vars_working$NBOT)
hist(all_vars_working$NANC)
hist(all_vars_working$NCUR)


#~~ Neutral model

folders <- list.files("neutral_scenario/fsc_run/", pattern = "fscrun[0-9]")

# collect likelihoods 
collect_vars_neutral <- function(folder){
  fsc_pars <- read_delim(paste0("neutral_scenario/fsc_run/", folder, "/afs/afs.bestlhoods"), delim = "\t")
}

all_vars_neutral <- folders %>% 
  map(safely(collect_vars_neutral)) %>% 
  purrr::transpose()

all_vars_working_neutral <- do.call(rbind, all_vars_neutral[[1]]) %>% 
  dplyr::mutate(lh_diff = MaxEstLhood - MaxObsLhood)

all_vars_working_neutral[which.max(all_vars_working_neutral$MaxEstLhood), ]

hist(all_vars_working_neutral$NCUR)
hist(all_vars_working_neutral$NANC)



#~~ AIC

# bot ML
all_vars_working[which.max(all_vars_working$MaxEstLhood), ]$MaxEstLhood

# neutral ML
all_vars_working_neutral[which.max(all_vars_working_neutral$MaxEstLhood), ]$MaxEstLhood


bot_AIC <- -2*(all_vars_working[which.max(all_vars_working$MaxEstLhood), ]$MaxEstLhood*log(10)) + 2*3
neutral_AIC = -2*(all_vars_working_neutral[which.max(all_vars_working_neutral$MaxEstLhood), ]$MaxEstLhood*log(10)) + 2*2

bot_AIC < neutral_AIC
