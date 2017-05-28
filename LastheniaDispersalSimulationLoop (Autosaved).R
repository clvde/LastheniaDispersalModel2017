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
Rmax_green <- 50 
Rmax_red <- 0
zmax_green <- 0
zmax_red <- 1
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

# Output a list of parameter values to file
parameters <- paste("meta_cols: ", meta_cols, ", ploidy: ", ploidy, ", disp_a_loci: ", disp_a_loci, ", disp_b_loci: ", disp_b_loci, ", env_loci: ", env_loci, ", neut_loci: ", neut_loci, ", Rmax_good: ", Rmax_good, ", Rmax_bad: ", Rmax_bad, ", nstar: ", nstar, ", p_mut: ", p_mut, ", sigma_mut: ", sigma_mut, ", nbhd_width: ", nbhd_width, ", env_length: ", env_length, ", t_max: ", t_max, ", init_loc_mean: ", init_loc_mean, ", k: ", k, ", disp_a_allele: ", disp_a_allele, ", disp_b_allele: ", disp_b_allele, ", env_allele: ", env_allele)
write(parameters, file = paste(filepathspec,"/parameters_from_sim.txt", sep = ''))

# ---------------------------------------------------------
# Simulation

# Create generation 0

current_population <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)

# Write it to file
write_name <- paste(filepathspec, 0, ".csv", sep="")
write.csv(current_population, write_name)

# Begin simulation loop
for (t in 1:t_max){
	print(paste("generation: ",t))
	
	# Make a dataframe to hold offspring. Budgets 10 times the "max carrying capacity" in a neighbourhood, nstar. 
	next_generation <- make_popn_dataframe(t, nstar*10, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	
	# Tracker which gives out unique IDs to individuals within a generation
	next_gen_ID_tracker <- 1
	
	#
	for (i in 1:nrow(current_population)){
		mom <- current_population[i,]
		Rmax <- environment(mom$location, Rmax_green, Rmax_red, t, env_length, env_change_speed)
		mates <- current_population[-i,]
		n_offspring <- reproduce(mom, nstar, Rmax, zmax, k, mates)
		dads_list <- matefinder1D(n_offspring, mom, mates, nbhd_width)
		dads_list_reformat <- convert_dads_list(dads_list)
		
		if (n_offspring > 0) {
			for (n in 1:n_offspring){
				dad <- dads_list_reformat[n]
				offspring <- make_offspring(mom, dad, t, next_gen_ID_tracker, p_mut)
				next_generation[next_gen_ID_tracker,] <- offspring
				next_gen_ID_tracker <- next_gen_ID_tracker + 1
			}
		}
	}

	# Get rid of empty rows in offpsring data frame, and then get the population size. Empty rows are identified by having a 0 in the ID column (column 2).
	next_generation <- next_generation[next_generation[,2]>0,]
	# Get population sizes
	print(paste("pop size: ",nrow(next_generation)))
	# This step creates discrete generations

	# write the parental generation to file before erasing them (annuals)
	write_name <- paste(filepathspec,t,".csv", sep="")
	write.csv(current_population, write_name)
	
	current_population <- next_generation
	
	if (nrow(current_population) == 0){
		print("extinction!")
		break
	}
}