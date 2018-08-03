# collect vars from bootstrap replicates

library(readr)
library(purrr)
library(stringr)
library(parallel)
library(dplyr)
library(skimr)
library(ggsci)
library(readr)
library(gridExtra)
library(ggthemr)
library(ggplot2)
source("sci_num.R")
options(scipen=999)



# how many bootstraps
nboot <- 200
num_sims <- 20
setwd("/home/emily/angsd_pilon")

# folders to collect from
# bootstrap folders
folders_boot <- list.files("nonpar_bootstrap/", pattern = "boot")

# simulations per bootstrap replicate folders
# folders_sim <- list.files("fastsimcoal_analyses/analyse_bootstrap_sfs/boot_1/", pattern = "fscrun[0-9]")
folders_sim <- paste0("fscrun", 1:num_sims)
# combine
folders_df <- data.frame("boot" = rep(folders_boot, each = length(folders_sim)), "sims" = rep(folders_sim, times = length(folders_boot)))

# collect likelihoods 
collect_vars <- function(folder_boot, folder_sim){
  fsc_pars <- read_delim(paste0("nonpar_bootstrap/", folder_boot, "/", folder_sim, "/afs/afs.bestlhoods"), delim = "\t")
}

# collect everything (usually just take maximum likelihood per bootstrap)
all_vars <- map2(folders_df$boot, folders_df$sims, possibly(collect_vars, NA_real_)) 
# put all in df
all_vars <- all_vars[!is.na(all_vars)]
all_vars_working <- bind_rows(all_vars)


system("mkdir data")
system("mkdir data/fsc_out")

write_delim(all_vars_working, path = "data/fsc_out/nonpar_boot_talk.txt", delim = " ")

hist(all_vars_working$NANC) # breaks == 100
hist(all_vars_working$NBOT) # breaks == 100
hist(all_vars_working$NCUR) # breaks == 100

par_ests <- all_vars_working %>% 
  mutate(boots = as.factor(MaxObsLhood)) %>% 
  group_by(boots) %>% 
 # arrange(desc(MaxEstLhood)) %>%  #   
  arrange(desc(MaxEstLhood), .by_group = TRUE) %>% 
  top_n(n = 2, MaxEstLhood) %>% # take top for each bootstrap
  mutate(NANC_dip = NANC/2,
         NBOT_dip = NBOT/2,
         NCUR_dip = NCUR/2)

hist(par_ests$NANC, breaks = 10)
hist(par_ests$NBOT, breaks = 10)
hist(par_ests$NCUR, breaks = 10)

#skimr::skim(all_vars_working)
#quantile(all_vars_working$NCUR)

#~~ 95% CIs

CI <- 0.95 # Check this

NANC_CI <- stats::quantile(par_ests$NANC_dip, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
NBOT_CI <- stats::quantile(par_ests$NBOT_dip, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
NCUR_CI <- stats::quantile(par_ests$NCUR_dip, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)


summary(par_ests)

#~~ Collect best empirical params

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
  dplyr::mutate(lh_diff = MaxEstLhood - MaxObsLhood,
                NANC_dip = NANC/2,
                NBOT_dip = NBOT/2,
                NCUR_dip = NCUR/2)

all_vars_working[which.max(all_vars_working$MaxEstLhood), ]

#~~ Plot non para boots

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")

mypal = pal_npg("nrc")(9)
mypal


# NANC
library(scales)

NANC_plot <- 
  ggplot(par_ests, aes(NANC_dip)) + 
  geom_histogram(colour = "grey45", fill = "grey45", binwidth = 100) +
  geom_vline(aes(xintercept = all_vars_working[which.max(all_vars_working$MaxEstLhood), ]$NANC_dip), 
             size = 1.5, linetype = "solid", show.legend = T, col = "#3C5488FF") +
  geom_errorbarh(aes(xmin = NANC_CI[1] , xmax = NANC_CI[2], y = 6),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  labs(y = "Frequency", x = expression(paste("Ancestral ", italic(N[e])))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  scale_x_continuous(labels = scales::unit_format("K", 1e-3, digits = 3)) +
  ggtitle('(A)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


# NBOT

NBOT_plot <-
  ggplot(par_ests, aes(NBOT_dip)) + 
  geom_histogram (fill = "grey45", binwidth = 25) +
  geom_vline(aes(xintercept = all_vars_working[which.max(all_vars_working$MaxEstLhood), ]$NBOT_dip), 
             size = 1.5, linetype = "solid", show.legend = T, col = "#2ca25f") +
  geom_errorbarh(aes(xmin = NBOT_CI[1] , xmax = NBOT_CI[2], y = 6),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  labs(y = "Frequency", x = expression(paste("Bottleneck ", italic(N[e])))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  ggtitle('(B)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


# NCUR

NCUR_plot <- 
  ggplot(par_ests, aes(NCUR_dip)) + 
  geom_histogram(fill = "grey45", bins = 15) +
  geom_vline(aes(xintercept = all_vars_working[which.max(all_vars_working$MaxEstLhood), ]$NCUR_dip), 
             size = 1.5, linetype = "solid", show.legend = T, col = "#E64B35FF") +
  geom_errorbarh(aes(xmin = NCUR_CI[1] , xmax = NCUR_CI[2], y = 6),
                 size = 0.8, color = "black", linetype = "solid", height = 0) +
  labs(y = "Frequency", x = expression(paste("Current ", italic(N[e])))) +
  theme(axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) +
  scale_x_continuous(labels = scales::unit_format("K", 1e-3, digits = 3)) +
  ggtitle('(C)') + theme(plot.title=element_text(hjust=0, size = 18, face = "plain"))


png("figs/Figure_3.png", units = "in", res = 300, width = 6, height = 15) # 22, 6

grid.arrange(NANC_plot, NBOT_plot, NCUR_plot,
             ncol = 1)

dev.off()

#~~ Summary stats

NANC_CI
NBOT_CI
NCUR_CI


all_vars_working[which.max(all_vars_working$MaxEstLhood), ]
