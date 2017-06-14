# This script contains parameter assignments and the simulation loop
# THIS SIMULATION PERFORMS A BURN IN

# File Path - where do you want simualtion files to go?
filepathspec = "/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Results/May 28th 2017 - Burn In/"

#Packages
library(statmod)

# Parameters - explained in LastheniaDispersalSimulationFunctions script
meta_cols <- 3 
meta_col_names <- c('generation','individual_ID','location')  
ploidy <- 2
disp_a_loci <- 5
disp_b_loci <- 5
env_loci <- 5
neut_loci <- 5
total_genome_length <- ploidy*(disp_a_loci+disp_b_loci+env_loci+neut_loci)
Rmax_green <- 50 
Rmax_red <- 50
zmax_green <- -1
zmax_red <- 1
nstar <- 20
p_mut <- 0.02
sigma_mut <- 0.0015 
nbhd_width <- 1 
env_length <- 1000 
t_max <- 100
env_change_speed <- 0
init_loc_mean <- 0
k <- 1
disp_a_allele <- 1 
disp_b_allele <- 2
env_allele <- -1

# Derived Params
disp_a_locus_1 <- meta_cols+1
disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
disp_b_locus_1 <- disp_a_locus_last + 1
disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
env_locus_1 <- disp_b_locus_last + 1
env_locus_last <- env_locus_1 + env_loci*ploidy - 1
neut_locus_1 <- env_locus_last + 1
neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

# Output a list of parameter values to file
# sink(paste(filepathspec,"/params.txt", sep = ''), type = "output")
# params <- as.data.frame(matrix(c(meta_cols, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci, total_genome_length, Rmax_green, Rmax_red, zmax_green, zmax_red, nstar, p_mut, sigma_mut, nbhd_width, env_length, t_max, env_change_speed, init_loc_mean, disp_a_allele, disp_b_allele, env_allele), nrow = 1, ncol = 22))
# colnames(params) <- c("meta_cols", "ploidy", "disp_a_loci", "disp_b_loci", "env_loci", "neut_loci", "total_genome_length", "Rmax_green", "Rmax_red", "zmax_green", "zmax_red", "nstar", "p_mut", "sigma_mut", "nbhd_width", "env_length", "t_max", "env_change_speed", "init_loc_mean", "disp_a_allele", "disp_b_allele", "env_allele")
# params
# sink()

# ---------------------------------------------------------
# Simulation

# Create generation 0

current_population <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)

#-----------------------------------------------------------
## try functions

makepopn <- make_popn_dataframe(0, nstar, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
makepopn
pop <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
pop
env <- environment(0, Rmax_green, Rmax_red, 0, env_length, env_change_speed)
env
Rmaxi <- env[1]
Rmaxi
zmaxi <- env[2]
zmaxi
indivi <- pop[1,]
indivi
zwi <- zw(indivi, env_locus_1, env_locus_last)
zwi
ldens <- localdensity(indivi, pop)
ldens 
ldens2 <- localdensity2(indivi, pop)
ldens2
Ri <- R(Rmaxi, zmaxi, k, zwi, ldens)
Ri
expoff <- expoffspring(zwi, ldens, 10000, Rmaxi, zmaxi, k)
expoff
noff <- reproduce(indivi, nstar, Rmaxi, zmaxi, k, pop)
noff
dads <- matefinder1D(noff, indivi, pop, nbhd_width)
dads
colnames(dads)
dadsvec <- convert_dads_list(dads)
dadsvec
dad <- pop[dadsvec[1],]
dad

offspring_list <- make_offspring(indivi, dad, 1, 1, p_mut)
offspring_list2 <- make_offspring2(indivi, dad, 1, 1, p_mut)
offspring_list
offspring_list2

zdai <- zda(indivi, disp_a_locus_1, disp_a_locus_last)
zdbi <- zdb(indivi, disp_b_locus_1, disp_b_locus_last)
offloc <- disperse1D(indivi$location, zdai, zdbi)
offloc








