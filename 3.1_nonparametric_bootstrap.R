# non-parametric bootstrap
# uncertainty in sites

library(readr)
library(purrr)
library(stringr)
library(parallel)
options(scipen=999)

# n sims for each np boot -- make 50
num_sim <- 20

setwd(paste0("/home/emily/angsd_pilon"))

# create non-parametric bootstrap folder
if (!dir.exists("/home/emily/angsd_pilon/nonpar_bootstrap")) {
  system("mkdir /home/emily/angsd_pilon/nonpar_bootstrap")
}

# create non-parametric bootstrap files folder
if (!dir.exists("/home/emily/angsd_pilon/nonpar_bootstrap_files")) {
  system("mkdir /home/emily/angsd_pilon/nonpar_bootstrap_files")
}

# read all bootstrapped sfs from angsd
all_sfs <- read_lines("afs_allsites_boot.sfs")

# how many bootstraps? 
nboot <- length(all_sfs)

# nonpar_bootstrap folder
nonpar_bootstrap_folder <- "/home/emily/angsd_pilon/nonpar_bootstrap_files/"

# delete all former runs
if (length(list.files(nonpar_bootstrap_folder)) != 0) { 
  system(paste0("rm -r ", nonpar_bootstrap_folder, "*"))
}

create_boot_sfs <- function(all_sfs_line, num_boot, folder) {
  # trim whitespace
  temp <- str_trim(all_sfs_line)
  # convert to numeric
  sfs <- as.numeric(unlist(str_split(temp, " ")))
  # add zeros for compatability with fsc
  sfs <- c(sfs, rep(0, length(sfs)-1))
  # create names
  sfs_names <- sapply(0:(length(sfs)-1), function(x) paste0("d0_", x))
  system(paste0("mkdir ", folder,  "boot_", num_boot))
  sink(paste0(folder, "boot_", num_boot, "/", "afs_MAFpop0.obs"))
  cat("1 observations")
  cat("\n")
  cat(sfs_names)
  cat("\n")
  cat(sfs)
  sink()
}

# create all nonparametric bootstrapped sfs
purrr::map2(all_sfs, 1:nboot, create_boot_sfs, nonpar_bootstrap_folder)

# nonpar_bootstrap folder
nonpar_bootstrap_run <- "/home/emily/angsd_pilon/nonpar_bootstrap/"

# delete all former runs
if (length(list.files(nonpar_bootstrap_run)) != 0) { 
  system(paste0("rm -r ", nonpar_bootstrap_run, "*"))
}

# run fsc
run_fsc <- function(run_num, num_boot){
  # create directory for fsc run
  system(paste0("mkdir fscrun", run_num))
  # paste all relevant files into directory
  system(paste0("cp ../../fastsimcoal_files/afs* fscrun", run_num))
  system(paste0("cp ../../nonpar_bootstrap_files/boot_", num_boot, "/afs_MAFpop0.obs fscrun", run_num))
  # change to directory
  setwd(paste0("/home/emily/angsd_pilon/nonpar_bootstrap/boot_", num_boot, "/fscrun", run_num))
  # run fsc
  system("/home/emily/programs/fsc26 -t afs.tpl -n 100000 -m -e afs.est -M -L 100 -q -w 0.01 -x --foldedSFS -C 5 --nosingleton")
  # change back
  setwd(paste0("/home/emily/angsd_pilon/nonpar_bootstrap/boot_", num_boot))
}

# run all
for (num_boot in 1:nboot) {
  
  if (!dir.exists(paste0("/home/emily/angsd_pilon/nonpar_bootstrap/boot_", num_boot))) {
    system(paste0("mkdir /home/emily/angsd_pilon/nonpar_bootstrap/boot_", num_boot))
  }
  setwd(paste0("/home/emily/angsd_pilon/nonpar_bootstrap/boot_", num_boot))
  
  cl <- makeCluster(getOption("cl.cores", 15)) # n cores
  parLapply(cl, 1:num_sim, run_fsc, num_boot)
  stopCluster(cl)
  
}

setwd("/home/emily/angsd_pilon/")