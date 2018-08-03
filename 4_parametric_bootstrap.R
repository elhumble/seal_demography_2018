# collect best lhoods from different fastsimcoal runs
# based on these parameters: simulate new sfs
# uncertainty in coalescent simulations
library(readr)
library(purrr)
library(stringr)
library(parallel)
library(dplyr)
library(ggplot2)
options(scipen=999)

setwd(paste0("/home/emily/angsd_pilon"))

# how many bootstraps
# nboot <- 100
# num_sim <- 10
# create SFS?
# create_SFS <- FALSE
# run bootstrap
# run_bootstrap <- TRUE

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



# determine highest likelihood
max_index <- which.max(all_vars_working$MaxEstLhood)
# min_index <- which.min(all_vars_working$lh_diff)
# folder with best lh 
folder_name <- folders[max_index]



### Parametric bootstrapping
sfs_sim_folder <- "/home/emily/angsd_pilon/sfs_simulation/"
  
# (1) manipulate par_file for data simulation
par_file <- readLines(paste0("fsc_run/", folder_name, "/afs/afs_maxL.par"))
# Number of sites in observed SFS: Number of lines in fasta file / 2 (roughly)
# 200000 DNA fragments of lenght 500
num_sites <- 36200
# new par_file
loci_row <- which(str_detect(par_file, "Number of independent loci")) + 1
datatype_row <- which(str_detect(par_file, "FREQ"))
# replace loci
par_file[loci_row] <- str_c(num_sites, " 0")
# replace datatype
par_file[datatype_row] <- "DNA 500 0  0.000000025 OUTEXP"
# write to new folder
if (!dir.exists(sfs_sim_folder)) {
  system(paste0("mkdir ", sfs_sim_folder)) 
}
writeLines(par_file, paste0(sfs_sim_folder, "afs.par"))
  
# (2) create 100 SFS in the bootstrap folder using fastsimcoal
nboot <- 100
setwd(sfs_sim_folder)
system(paste("/home/emily/programs/fsc26 -i afs.par -n ", nboot, " -j -m -s0 -x -q", sep = ""))
setwd(paste0("/home/emily/angsd_pilon"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Compare SFS distributions     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Import empirical sfs

emp_sfs <- readLines("fastsimcoal_files/afs_MAFpop0.obs")[3]
emp_sfs <- as.numeric(unlist(str_split(emp_sfs, " ")))%>%
   .[1:71]
  
  
#~~ Function to import simulated sfs

import_sfs <- function(file) {
    sfs <- readLines(file)[3]
    sfs <- as.numeric(unlist(str_split(sfs, "\t"))) %>%
      .[1:71]
    }
  
  
#~~ Read in all simulated sfs

nsim <- 100
sim_sfs <- list()

for (i in 1:nsim) {
  sim_sfs[i] <- lapply(paste0("sfs_simulation/afs/afs_", i,"/afs_MAFpop0.obs"), import_sfs)
}
  
sim_sfs <- as.data.frame(matrix(unlist(sim_sfs), nrow=length(unlist(sim_sfs[1]))))
  
# CIs
CI <- 0.95
  
CIs <- apply(sim_sfs, 1, function(x) stats::quantile(x, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE))
CI5 <- CIs[1,]
CI95 <- CIs[2,]
sim_sfs_means <- rowMeans(sim_sfs)



#~~ Create dataframe for plotting

sfs <- data.frame(emp_sfs) %>%
  mutate(CI5 = CI5,
         CI95 = CI95,
         sim_mean = sim_sfs_means,
         sfs = seq(1:71)) %>%
  cbind(sim_sfs)



#~~ Create plotting theme

library(ggthemr)

ggthemr(palette = "pale", layout = "clean", 
        line_weight = 0.7, text_size = 20, type = "outer")



library(ggsci)
mypal = pal_npg("nrc")(9)
mypal

a <- ggplot(sfs) +
  geom_col(aes(x = sfs, y = emp_sfs), fill = "#3C5488FF", alpha = 0.8) +
  geom_point(aes(x = sfs, y = sim_mean), size = 0.1, col = "grey30") +
  geom_linerange(aes(ymin = CI5, ymax = CI95, x = sfs), col = "grey30") +
  geom_col(aes(x = sfs, y = sim_mean), fill = "grey30", alpha = 0.5) +
  coord_cartesian(ylim=c(0,16000)) +
  labs(y = "Number of SNPs", x = "Allele frequency") +
  scale_y_continuous(labels = scales::unit_format("K", 1e-3, digits = 1)) +
  theme(plot.margin=unit(c(1,1,1,1),"lines"),
        axis.text.x = element_text(face = "plain"),
        axis.text.y = element_text(face = "plain"),
        axis.title.x = element_text(face = "plain"),
        axis.title.y = element_text(face = "plain")) 

# monomorphic
b <- ggplot(sfs) +
  geom_col(aes(x = sfs, y = emp_sfs), fill = "#3C5488FF", alpha = 0.8) +
  geom_point(aes(x = sfs, y = sim_mean), size = 0.1, col = "grey30") +
  geom_linerange(aes(ymin = CI5, ymax = CI95, x = sfs), col = "grey30") +
  coord_cartesian(ylim=c(15000000,20000000)) +
  scale_y_continuous("", labels=scales::unit_format("M", 1e-5, digits = 3), breaks = c(16000000,20000000)) +
  geom_col(aes(x = sfs, y = V1), fill = "grey30", alpha = 0.5) +
  #theme(plot.margin=unit(c(0,0.5,-0.5,0), "lines"),
        theme(plot.margin=unit(c(1,1,-0.5,1), "lines"),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")

# doubletons
c <- ggplot(sfs) +
  geom_col(aes(x = sfs, y = emp_sfs), fill = "#3C5488FF", alpha = 0.8) +
  geom_point(aes(x = sfs, y = sim_mean), size = 0.1, col = "grey30") +
  geom_linerange(aes(ymin = CI5, ymax = CI95, x = sfs), col = "grey30") +
  coord_cartesian(ylim=c(65000,67000)) +
  scale_y_continuous("", labels=scales::unit_format("K", 1e-3, digits = 3), breaks = c(66000)) +
  geom_col(aes(x = sfs, y = V1), fill = "grey30", alpha = 0.5) +
  #theme(plot.margin=unit(c(0,0.5,-0.5,0), "lines"),
  theme(plot.margin=unit(c(1,1,-0.5,1), "lines"),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none")
  
library(gridExtra)


gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
gC <- ggplotGrob(c)

maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gC$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)

grid.arrange(gB,gC, gA, ncol=1, heights=c(0.15,0.1,0.85))


png("/home/emily/angsd_pilon/figs/Figure_2.png", units = "in", res = 300, width = 8, height = 7)
grid.arrange(gB,gC, gA, ncol=1, heights=c(0.15,0.1,0.85))
dev.off()


 
# #### Second part
# 
# setwd(paste0("/home/martin/nes/nes_RAD"))
# # (3) Run Fsc 50 times in each folder
# # create new folder for likelihood maximization of all 100 bootstrap SFS samples
# if (!dir.exists("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")) {
#   system("mkdir /home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")
# }
# setwd("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")
# 
# if (length(list.files("./")) != 0) { 
#   system("rm -r *")
# }
# 
# run_fsc <- function(run_num, num_boot){
#   # create directory for fsc run
#   system(paste0("mkdir fscrun", run_num))
#   # paste all relevant files into directory
#   system(paste0("cp ../../fastsimcoal_files/nes* fscrun", run_num))
#   system(paste0("cp ../../bootstrap/nes/nes_", num_boot, "/nes_MAFpop0.obs fscrun", run_num))
#   # change to directory
#   setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot, "/fscrun", run_num))
#   # run fsc
#   system("~/bin/fsc26 -t nes.tpl -n 10000 -m -e nes.est -M -L 40 -q -w 0.01 -x --foldedSFS -C 5 --nosingleton")
#   # change back
#   setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
# }
# 
# # run all
# for (num_boot in 1:nboot) {
#   
#   if (!dir.exists(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))) {
#     system(paste0("mkdir /home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
#   }
#   setwd(paste0("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs/boot_", num_boot))
#   
#   cl <- makeCluster(getOption("cl.cores", 25))
#   parLapply(cl, 1:num_sim, run_fsc, num_boot)
#   stopCluster(cl)
#   
#   setwd("/home/martin/nes/nes_RAD/fastsimcoal_analyses/analyse_bootstrap_sfs")
# }
# 
# setwd("/home/martin/nes/nes_RAD/")
# 
# 
# 
# 
# 
