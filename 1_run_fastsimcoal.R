# run several iterations of fastsimcoal

library(parallel)

if (!dir.exists("/home/emily/angsd_pilon/fsc_run")) {
  system("mkdir /home/emily/angsd_pilon/fsc_run")
}

setwd("/home/emily/angsd_pilon/fsc_run")
if (length(list.files("./")) != 0) { 
  system("rm -r *")
}

# system("~/bin/fsc26")
#folders <- list.files("fsc_run/", pattern = "fsc_run[0-9]")

run_fsc <- function(run_num){
  # create directory for fsc run
  system(paste0("mkdir fscrun", run_num))
  
  # paste all relevant files into directory
  system(paste0("cp ../fastsimcoal_files/afs* fscrun", run_num))
  
  # change to directory
  setwd(paste0("/home/emily/angsd_pilon/fsc_run/fscrun", run_num))
  
  #if (run_num < 51) {
  # write fsc_run file 
  writeLines(c(paste0("/home/emily/angsd_pilon/fsc_run/fscrun", run_num),
               "-t afs.tpl -n 100000 -m -e afs.est -M -L 100 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton"),
             paste0("/home/emily/angsd_pilon/fsc_run/fscrun", run_num, "/fsc_run.txt"))
  
  #}
  
  # run fsc
  system("~/programs/fsc26")
  
  # change back
  setwd(paste0("/home/emily/angsd_pilon/fsc_run/"))
  
}

# run all
cl <- makeCluster(getOption("cl.cores", 15))
parLapply(cl, 1:100, run_fsc) # number of simulations (increase)
stopCluster(cl)

setwd("/home/emily/angsd_pilon/")

