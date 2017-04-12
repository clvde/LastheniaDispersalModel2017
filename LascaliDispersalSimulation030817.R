# Beginning of Dispersal Simulations
# Author: Courtney Van Den Elzen
# Rotation 3 Project: Melbourne and Flaxman Labs
# February 2017 - May 2017

# The goal of this project is to set up a simulation modelling framework which can be used to describe the population dynamics and spread of annual, wind-dispersed grassland plants

# MAP OF THIS SCRIPT:

# (1) Assumptions
# (2) Constants
# (3) Population Data Frame Creation
# (4) Reproductive Functions
# (5) Dispersal Functions
# (6) The Environment
 
# -------------------------------- (1) Assumptions: NEED TO COMPLETE THIS SECTION -------------------------------- 
	# Discrete generations
	# self-incompatible
	# diploid organisms
	# 5 loci controlling dispersal evolution (there are multiple traits that contribute to the dispersal phenotype)
	# 5 loci controlling environmental evolution (there are multiple traits that contribute to the dispersal phenotype)
	# 5 neutral loci as controls - allow us to spatial sorting and change in neutral genetic diversity over space etc.
	# no new mutations
	# all seeds from a given mom are fathered by the same dad (this is probably unrealistic but we can change it after)
	# hermaphroditic individuals
	# free recombination (every locus is on a different chromosome). Probability of recombination constant between every locus

   
# -------------------------------- (2) CONSTANTS -------------------------------- 
# This should be changed into a simple list of the constants and their definitions. The assignments of values to these constants should only be below. 

# meta_cols: The number of metadata columns present. In the first iteration of this simulation, there are 3: (1) Generation (2) Individual ID (a unique (within generation) identifier of the individual) (3) Location
# ploidy: The ploidy level of the genome. In general this will always be equal to 2 (diploid organisms)
# disp_a_loci: The number of diploid loci that contribute to the dispersal parameter a
# disp_b_loci: The number of diploid loci that contribute to the dispersal parameter b
# env_loci: The number of diploid loci that contribute to the environmental (or fitness) parameter.
# neut_loci: The number of diploid loci that are neutral in the genome and do not contribute to any phenotype 
# total_genome_length: The total number of loci in the genome, but one diploid locus contributes 2. i.e. t_g_l <- ploidy*(disp_a_loci+disp_b_loci+env_loci+neut_loci)
# Rmax_good: The intrinsic growth rate in the "good" habitat/environment 
# Rmax_bad: The intrinsic growth rate in the "bad" habitat/environment
# nstar: The maximum carrying capacity, or the total reproductive output of a neighbourhood
# p_mut: The probability of acquiring a mutation. We assume that only one mutation max is acquired per individual because the porbability of more than one is sufficiently sml
# sigma_mut: The standard deviation of the mutation kernel - mutations are drawn from a normal distribution with mean p_mut and sd sigma_mut
# nbhd_width: The "width" or "size" or "standard deviation" of the neighbourhood. The neighbourhood weight is a normal distribution with mean of the maternal location & sd nbhd_width
# env_length: The length of the "good" habitat location. This should be varied (i.e. the neighbourhood size to environment size varies) as this may have consequences for dipsersal evo

# -------------------------------- (3) POPULATION DATA FRAME CREATION -------------------------------- 

# Structure of the data frame which holds all of the population information (number of metadata columns may change):
   # Column 1: generation
   # Column 2: individual ID (number from 1 to pop size in the current generation)
   # Column 3: location
   # Column 4-13: dispersal trait 1 (a = mean) loci (2 columns per locus because of diploidy)
   # Column 14-23: dispersal trait 2 (b = shape) loci 
   # Column 24-33: environmental loci
   # Column 34-43: neutral loci
   
# Makes the original data frame with the right numbers of columns and the right column labels -- use this to create an empty data frame to store the data from each generation

make_popn_dataframe <- function(t, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci){
	current_population <- data.frame(x=0)
	current_population[,1:meta_cols] <- c(0)
	colnames(current_population) <- meta_col_names

# 	disp_a_locus_1 <- meta_cols+1
	disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
	disp_b_locus_1 <- disp_a_locus_last + 1
	disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
	env_locus_1 <- disp_b_locus_last + 1
	env_locus_last <- env_locus_1 + env_loci*ploidy - 1
	neut_locus_1 <- env_locus_last + 1
	neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

	# creates all of the necessary new columns for the loci
	for (i in seq(from=meta_cols, to=(total_genome_length+meta_cols), by=1)){
		current_population[,i] <- 0
	}

	# Names the columns properly

	for (i in (meta_cols+1):ncol(current_population)){
	
		if (meta_cols%%2 != 0){
	
			if (i%%2 == 0) {
				if (4 <= i && i <= disp_a_locus_last){
					j = i-meta_cols
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispa_locus',k,'1',sep = "_")
				} else if (disp_b_locus_1 <= i && i <= disp_b_locus_last){
					j = i-disp_a_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispb_locus',k,'1',sep = "_")
				} else if (env_locus_1 <= i && i <= env_locus_last) {
					j = i-disp_b_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('env_locus',k,'1',sep = "_")
				} else if (neut_locus_1 <= i && i <= neut_locus_last) {
					j = i-env_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('neut_locus',k,'1',sep = "_")
				}
			} else if (i%%2 != 0){
				if (4 <= i && i <= disp_a_locus_last){
					j = i-meta_cols
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispa_locus',k,'2',sep = "_")
				} else if (disp_b_locus_1 <= i && i <= disp_b_locus_last){
					j = i-disp_a_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispb_locus',k,'2',sep = "_")
				} else if (env_locus_1 <= i && i <= env_locus_last) {
					j = i-disp_b_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('env_locus',k,'2',sep = "_")
				} else if (neut_locus_1 <= i && i <= neut_locus_last) {
					j = i-env_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('neut_locus',k,'2',sep = "_")
				}
			}
		} else if (meta_cols%%2 == 0) {
		
			if (i%%2 != 0) {
				if (4 <= i && i <= disp_a_locus_last){
					j = i-meta_cols
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispa_locus',k,'1',sep = "_")
				} else if (disp_b_locus_1 <= i && i <= disp_b_locus_last){
					j = i-disp_a_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispb_locus',k,'1',sep = "_")
				} else if (env_locus_1 <= i && i <= env_locus_last) {
					j = i-disp_b_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('env_locus',k,'1',sep = "_")
				} else if (neut_locus_1 <= i && i <= neut_locus_last) {
					j = i-env_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('neut_locus',k,'1',sep = "_")
				}
			} else if (i%%2 == 0){
				if (4 <= i && i <= disp_a_locus_last){
					j = i-meta_cols
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispa_locus',k,'2',sep = "_")
				} else if (disp_b_locus_1 <= i && i <= disp_b_locus_last){
					j = i-disp_a_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('dispb_locus',k,'2',sep = "_")
				} else if (env_locus_1 <= i && i <= env_locus_last) {
					j = i-disp_b_locus_last
					k = round((j/2)+0.0000001)
					colnames(current_population)[i] <- paste('env_locus',k,'2',sep = "_")
				} else if (neut_locus_1 <= i && i <= neut_locus_last) {
					j = i-env_locus_last
					k = round((j/2)+0.0000001)
				colnames(current_population)[i] <- paste('neut_locus',k,'2',sep = "_")
				}	
			}
		}
	}
	current_population <- current_population[-1,]
	return(current_population)
}


# This uses the make_popn_dataframe function to create a data frame and then fill it with n identical individuals -- which is how the simulation starts. 
# Possible modifications to make: non-identical individuals, individuals at different locations

make_pop <- function(t, n, init_location, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci){
	
	curr_pop <- make_popn_dataframe(t, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	
	disp_a_genome <- as.vector(matrix(disp_a_allele/(ploidy*disp_a_loci), nrow = 1, ncol = ploidy*disp_a_loci))
	disp_b_genome <- as.vector(matrix(disp_b_allele/(ploidy*disp_b_loci), nrow = 1, ncol = ploidy*disp_b_loci))
	env_genome <- as.vector(matrix(env_allele/(ploidy*env_loci), nrow = 1, ncol = ploidy*env_loci))
	neut_genome <- as.vector(matrix(0, nrow = 1, ncol = ploidy*neut_loci))
	
		
	for (i in 1:n){
		init_loc_rand <- rnorm(1, mean = init_location, sd = nbhd_width) # the initial location is random but all individuals will fall within the same neighbourhood (i.e. the standard deviation is the same as the neighbourhood width)

		curr_pop[i,] <- c(0, i, init_loc_rand, disp_a_genome, disp_b_genome, env_genome, neut_genome)
	}
	return(curr_pop)
}

# -------------------------------- (4) REPRODUCTIVE FUNCTIONS -------------------------------- 

# Calculates the fitness phenotype for an individual

zw <- function(individual, start_locus, end_locus){
	zfit <- sum(individual[start_locus:end_locus])
	zfit <- as.numeric(zfit)
	return(zfit)
}


# Calculates the expected fecundity of an individual

R <- function(Rmax, k, zw){
	Ri <- Rmax*exp(-k*(zw^2))
	return(Ri)
}


# Calculates the expected number of offspring (as in Phillips 2015)

expoffspring <- function(zwi, nix, nstar, Rmax, k){
	Ri <- R(Rmax,k,zwi)
	alpha <- (Rmax - 1)/nstar
	EO <- Ri/(1+alpha*nix)
	return(EO)
}

# A function to choose the dads of the offspring. This function does not determine the number of offspring, but instead accepts that as a parameter. My 

matefinder1D <- function(nbabies, mom, current_population, nbhd_width){
	# construct a probability density function (multinomial) of mating with every dad (by distance) and then sample from that based on mom's density-independent fecundtiy
	
	momID <- mom$individual_ID
	mom_loc <- mom$location
	dads <- current_population[-momID,] # keep in mind that taking mom out will make the row number different from the individual ID number
	dadIDs <- as.data.frame(dads$individual_ID)
	dadIDs[,2] <- dads$location
	dadIDs[,3] <- 0	
	colnames(dadIDs) <- c('individual_ID','location','probability')
	
	# fills in the probability of mating with each dad given the distance away, the 
	for (i in 1:nrow(dads)){
		prob <- dnorm(dadIDs$location[i], mom_loc, nbhd_width)
		dadIDs[i,3] <- prob
	}
	
	# now we have a list of probabilities for each dad, which parameterizes a multinomial distribution. 
	babies_vec <- rmultinom(1,size = nbabies, prob=dadIDs[,3])
	
	dadIDs[,4] <- babies_vec
	colnames(dadIDs)[4] <- "offspring_counts"
	
	return(dadIDs)	
}

convert_dads_list <- function(dadIDs){

	dad_indices <- which(dadIDs$offspring_counts>0, arr.ind=TRUE)[,1]
	IDs <- as.vector(dadIDs$individual_ID)[dad_indices]
	offspring <- as.vector(dadIDs$offspring_counts)[dad_indices]
	dads_by_offnum <- c()
	
	for (i in 1:length(IDs)){
		tempvec <- matrix(data = IDs[i], nrow = 1, ncol = offspring[i])
		dads_by_offnum <- c(dads_by_offnum, tempvec)
	}
	
	return(dads_by_offnum)
}

# A function to produce one child to a given mom. This is a helper function for reproduce ()

make_offspring <- function(mom, dad, generation, indiv_num){ # this whole function assumes that the first two columns of any indidividual have their generation number and their individual number. So locus 1 is in the third entry of the vector that defines every individual, and so on.
	
	# give the baby it's generation number and individual number
	baby_vec_length = meta_cols + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1] <- generation
	baby[,2] <- indiv_num
	baby[,3] <- disperse1D(mom[,3], zda(mom, disp_a_locus_1, disp_a_locus_last), zdb(mom, disp_b_locus_1, disp_b_locus_last))
	
	# create the baby's genome
	locus_vec <- seq(from = 3, to = total_genome_length+2, by = 2)
	for (i in locus_vec){
		if (i%%2 != 0){ # if this is an the first chromosome of a pair, inherit from mom at random
			momallele = round(runif(1)) + i # choose either the allele at locus i_1 or at locus i_2
			baby[,i] = mom[,momallele] 
		}	else { # if this is an odd chromosome, inherit from dad at random
				dadallele = round(runif(1))+i	# choose either the allele at locus i_1 or at locus i_2 from dad's genome
				baby[,i+1] = dad[,dadallele]	
			}
	}
	
	# now mutate if needed
	
	if (runif(1) < 1){
		focal_locus <- sample(c(4:total_genome_length+3),1) # sample a random locus
		mut_allele <- rnorm(1, baby[focal_locus], sigma_mut)
		baby[focal_locus] <- mut_allele
	}
	return(baby)
}

# every individual contributes to the local density based on their distance only -- so here we assume that the fitness phenotype does NOT effect the density that individual contributes. In other words, the every individual's effective density contribution is the same, regardless of its phenotype. 
# This also includes the focal individual in the calculations (this shouldn't matter as long as it's consistent wit nstar)

localdensity <- function(mom, current_population){
	
	distances <- abs(current_population[,3]-mom[3])
	local_density <- 0 
	
	for (i in 1:nrow(current_population)){
		focal_indiv <- current_population[i,]
		dix <- abs(mom$location - focal_indiv$location)
		focal_indiv_contribution <- dnorm(dix,mom$location,nbhd_width)
		local_density <- local_density + focal_indiv_contribution
	}
	return(local_density)
}

# Creates the actual number of offspring for a mother from its expected value (Poisson distributed around a mean of the expected)
reproduce <- function(mom, nstar, Rmax, k, current_population){
	# nstar is the population size that would be achieved if all individuals had Rmax fecundity. 
	# Nix is the current local population density around mom
	
	nix <- localdensity(mom, current_population)
	zwi <- as.numeric(zw(current_population, env_locus_1, env_locus_last))
	exp_offspring <- expoffspring(zwi, nix, nstar, Rmax, k)
	num_offspring <- rpois(1, exp_offspring)
	
	return(num_offspring)	
}


# -------------------------------- (5) DISPERSAL FUNCTIONS --------------------------------

# Calculates the first dispersal phenotype (a- mean distance) for an individual.

zda <- function(individual, start_locus, end_locus){
	# this assumes that you have one vector describing an individual. entries 13-22 describe the dispersal alleles
	zdispa <- sum(individual[start_locus:end_locus])
	return(zdispa)
}


zdb <- function(individual, start_locus, end_locus){
	# this assumes that you have one vector describing an individual. entries 13-22 describe the dispersal alleles
	zdispb <- sum(individual[start_locus:end_locus])
	return(zdispb)
}


disperse1D <- function(mom_loc, zda, zdb){
	# assumes mom's location is in one dimension
	dist <- rinvgauss(1, zda, zdb)
	dir <- sample(c(-1,1), 1)
	new_loc <- mom_loc + dir*dist
	return(new_loc)	
}


# A function to choose how far away a given seed disperses. a is the mean parameter, which I think I will define as the maternal dispersal phenotype. b controls the shape of the distribution. For now this will remain constant over time and not evolve... perhaps in later iterations of this model we could let this evolve freely, or not even impose a certain kernel! Maybe the distribution shoudl be an emergent property?! 

# this assumes 2D dispersal - with every individual having an x and a y coordinate
# could also just import mom's row and extract the phenotypes within the function

# disperse2D <- function(mom_loc, zda, zdb){
	
	# xm <- momloc[1] # assumes moms location
	# ym <- momloc[2]
	# theta <- 360*runif(1)
	# dist <- rinvgauss(1, zda, zdb)
	
	# if (theta == 0 || theta == 360) {
		
		# xb <- xm
		# yb <- ym + dist
		
	# } else if (0 < theta < 90) {
		
		# xb <- xm - dist*sin(theta)
		# yb <- ym + dist*cos(theta)
		
	# } else if (theta == 90) {
		
		# xb <- xm - dist
		# yb <- ym
		
	# } else if (90 < theta < 180) {
		
		# xb <- xm - dist*sin(180-theta)
		# yb <- ym - dist*cos(180-theta)
		
	# } else if (theta == 180) {
		
		# xb <- xm 
		# yb <- ym - dist
				
	# } else if (180 < theta < 270) {
		
		# xb <- xm + dist*sin(270-theta)
		# yb <- ym - dist*cos(270-theta)
		
	# } else if (theta == 270) {
		
		# xb <- xm + dist
		# yb <- ym
		
	# } else if (270 < theta < 360) {
		
		# xb <- xm + dist*sin(360-theta)
		# yb <- ym + dist*cos(360-theta)
		
	# }	
	
	# baby_loc <- c(xb,yb)
	
	# return(baby_loc)
	
# }


# -------------------------------- (6) THE ENVIRONMENT -------------------------------- 

# this function provides a map of what the Rmax is at any given point in space asa function of time
# we use the value of t as the middle of the good environment. So at t=0, the good environment extends from [-env_length/2, env_length/2].

environment <- function(dix, Rmax_good, Rmax_bad, t, env_length){
	
	low_lim <- t-(env_length/2)
	upp_lim <- t+(env_length/2)
	
	if (low_lim < dix && upp_lim > dix){
		return(Rmax_good)
	} else {
		return(Rmax_bad)
	}
}

# -------------------------------- (7) RUN SIMULATION LOOP --------------------------------

# Constants: assigning values

meta_cols <- 3 # the number of metadata columns in the matrix (generation number, individual number, location)
meta_col_names <- c('generation','individual_ID','location')
ploidy <- 2
disp_a_loci <- 5
disp_b_loci <- 5
env_loci <- 5
neut_loci <- 5
total_genome_length <- ploidy*(disp_a_loci+disp_b_loci+env_loci+neut_loci)
Rmax_good <- 50 # (50 seeds per inflorescence, 1 inflorescence MAX)
Rmax_bad <- 0
sigma <- 1 # effective size of spatial neighbourhood
nstar <- 1000
b <- 1 # shape parameter of the inverse gaussian distribution
p_recomb <- 0.0001
p_mut <- 0.00001 # probability of mutation at every allele
sigma_mut <- 0.001 # standard deviation of the mutation kernel
nbhd_width <- 1 # can set this equal to 1 without loss of generality as the whole environment can be large or small relative to this. 
env_length <- 10 # this should be varied so the gene flow (i.e. the neighbourhood size to environment size varies) as this may have consequences for dipsersal evo
t_max <- 10


# Make the initial population (generation 0) - nstar identical individuals all at location 0

current_population <- make_init_pop(nstar, init_location, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)

for (t in 1:t_max){
	# (1) Reproduction
	# (2) Parental Death
	# (3) Dispersal (but this is density independent)
	# (3) F1 Reproduction
	
	
	# (1) & (2) - offspring dispersal is built into the make_offspring function. 
	
	next_generation <- make_popn_dataframe(t,meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	next_gen_ID_tracker <- 1
	
	for (i in 1:nrow(current_population)){
		mom <- current_population[i,]
		Rmax <- environment(mom$location, Rmax_good, Rmax_bad, t, env_length)
		mates <- current_population[-mom,]
		n_offspring <- reproduce(mom, nstar, Rmax, k, mates)
		dads_list <- matefinder1D(n_offspring, mom, mates, nbhd_width)
		dads_list_reformat <- convert_dads_list(dads_list)
		
		for (n in 1:n_offspring){
			dad <- dads_list_reformat[n]
			new_babe <- make_offspring(mom, dad, current_population, t, indiv_num)
			next_generation[next_gen_ID_tracker,] <- new_babe
			next_gen_ID_tracker <- next_gen_ID_tracker + 1
		}
		
	}
	
	# write the parental generation to file before erasing them (annuals)
	write_name <- paste("/Users/Courtney/Documents/Simulation Practice Files/lascali_sim_dispevo_only_gen_",t,".csv"sep="")
	write.csv(current_population,write_name,col_names = TRUE, row_names = TRUE)
	
	current_population <- next_generation
	
}
	







