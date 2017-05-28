# This script contains parameter assignments and the simulation loop

# File Path - where do you want simualtion files to go?
filepathspec = "/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Practice Files"

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
Rmax_1 <- 50 
Rmax_2 <- 50
zmax_1 <- -1
zmax_2 <- 1
nstar <- 100
p_mut <- 0.00001 
sigma_mut <- 0.001 
nbhd_width <- 1 
env_length <- 10 
t_max <- 100
env_change_speed <- 0.1
init_loc_mean <- 0
k <- 1
disp_a_allele <- 1 
disp_b_allele <- 2
env_allele <- 0

# Derived Params
disp_a_locus_1 <- meta_cols+1
disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
disp_b_locus_1 <- disp_a_locus_last + 1
disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
env_locus_1 <- disp_b_locus_last + 1
env_locus_last <- env_locus_1 + env_loci*ploidy - 1
neut_locus_1 <- env_locus_last + 1
neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

#---------------------------------------------------------------
# Testing functions

popdf <- make_popn_dataframe(t_max, nstar, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
pop <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
indiv_i <- pop[1,]
current_population <- pop[-1,]
zwi <- zw(indiv_i, env_locus_1, env_locus_last)
Ri <- R(Rmax_1, zmax_1, k, -1, 1)
e0 <- 10 #expoffspring(zmax_1, 10, 2*nstar, Rmax_1, zmax_1, k)
noffspring <- 15 #reproduce(indiv_i, nstar, Rmax_1, zmax_1, k, current_population)
dadIDs <- matefinder1D(noffspring, indiv_i, current_population, nbhd_width)
dads_by_offnum <- convert_dads_list(dadIDs)
baby <- make_offspring(indiv_i, pop[dads_by_offnum[1],], 1, 1)
loc_dens <- localdensity(indiv_i, pop)
loc_dens 
indiv_i[,3]
indiv_i[3]
current_population[,3]-indiv_i[,3]



