# Check the model is a good fit to the data
options(scipen = 999)
library(stringr)
library(dplyr)

emp_sfs <- readLines("empirical_sfs/afs_MAFpop0.obs")[3]
emp_sfs <- as.numeric(unlist(str_split(emp_sfs, " ")))%>%
  .[1:71]

sim_sfs <- readLines("sfs_sim/afs_MAFpop0.obs")[3]
sim_sfs <- as.numeric(unlist(str_split(sim_sfs, "\t"))) %>%
  .[1:71]

sim_sfs_orig <- readLines("~/Desktop/afs_MAFpop0.obs")[3]
sim_sfs_orig <- as.numeric(unlist(str_split(sim_sfs_orig, "\t"))) %>%
  .[1:71]

barplot(emp_sfs[-c(1,2)])
barplot(sim_sfs[-c(1,2)])
barplot(sim_sfs_orig[-c(1,2)])

df <- data.frame(emp_sfs, sim_sfs, sim_sfs_orig) %>%
  mutate(sfs = seq(1:nrow(.)))


ggplot(df[-c(1,2),]) +
  geom_col(aes(x = sfs, y = emp_sfs), fill = "red", alpha = 0.5) +
  geom_col(aes(x = sfs, y = sim_sfs), fill = "blue", alpha = 0.5) +
  geom_col(aes(x = sfs, y = sim_sfs_orig), fill = "grey", alpha = 0.5)

# include monomorphic sites as axis cutoff

