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
Rmax_good <- 50 
Rmax_bad <- 0
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

current_population <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
current_population
curr_pop <- make_popn_dataframe(t, 1, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
curr_pop
disp_a_genome <- rep("a", ploidy*disp_a_loci)
disp_a_genome
disp_b_genome <- rep("b", ploidy*disp_b_loci)
disp_b_genome
env_genome <- rep("e", ploidy*env_loci)
env_genome
neut_genome <- rep("n", ploidy*neut_loci)
neut_genome

#curr_pop[c(1,2),] <- c("gen", "ID", "loc", disp_a_genome, disp_b_genome, env_genome, neut_genome)	

curr_pop[c(1:10),] <- 0
curr_pop[,1] <- "gen"
curr_pop[,2] <- "ID"
curr_pop[,3] <- rnorm(10, mean = 0, sd = 1)

curr_pop