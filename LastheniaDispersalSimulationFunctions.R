# Dispersal Evolution Simulations
# Author: Courtney Van Den Elzen
# Rotation 3 Project: Melbourne and Flaxman Labs
# February 2017 - June 2017

# The goal of this project is to set up a simulation modelling framework which can be used to describe the population dynamics and spread of annual, wind-dispersed grassland plants

# MAP OF THIS SCRIPT:

# (1) Assumptions
# (2) Constants
# (3) Population Data Frame Creation
# (4) Reproductive Functions
# (5) Dispersal Functions
# (6) The Environment
 
# -------------------------------- (1) ASSUMPTIONS -------------------------------- 

 # (1) Discrete generations
 # (2) Diploid organisms 
 # (3) Self-incompatible
 # (4) Two traits controlling dispersal evolution, both have multiple loci controlling them. First controls average dispersal distance, second controls shape param of kernel 
 # (5) One trait controlling environmental evolution. Has multiple loci controlling it.
 # (6) Multiple neutral loci used to observe neutral accumulation of genetic variation, as well as allele surfing
 # (7) Mutations arise with probability 0.00001 in any individual. It is assumed that 2 or more mutations per individual never happens, since probability is max 0.00001^2. 
 # (8) Multiple fathers possible - probability of paternity based entirely on distance from mother (probably roughly true for passively dispersed pollen)
 # (9) Every individual has chance at being both mother and father
 # (10) Free recombination (every locus is on a different chromosome). Implemented as randomly drawing one allele from each parent at each locus
 # (11) XXXXXXX - need to complete this section
   
# -------------------------------- (2) CONSTANTS -------------------------------- 

# meta_cols: The number of metadata columns present. In the first iteration of this simulation, there are 3: (1) Generation (2) Individual ID (a unique (within generation) identifier of the individual) (3) Location
# ploidy: The ploidy level of the genome. In general this will always be equal to 2 (diploid organisms)
# disp_a_loci: The number of diploid loci that contribute to the dispersal parameter a
# disp_b_loci: The number of diploid loci that contribute to the dispersal parameter b
# env_loci: The number of diploid loci that contribute to the environmental (or fitness) parameter.
# neut_loci: The number of diploid loci that are neutral in the genome and do not contribute to any phenotype 
# total_genome_length: The total number of loci in the genome, but one diploid locus contributes 2. i.e. t_g_l <- ploidy*(disp_a_loci+disp_b_loci+env_loci+neut_loci)
# Rmax_good: The intrinsic growth rate in the "good" habitat/environment 
# Rmax_bad: The intrinsic growth rate in the "bad" habitat/environment
# nstar: The interdensity an individual feels at which the expected number of 
# p_mut: The probability of acquiring a mutation. We assume that only one mutation max is acquired per individual because the porbability of more than one is sufficiently sml
# sigma_mut: The standard deviation of the mutation kernel - mutations are drawn from a normal distribution with mean p_mut and sd sigma_mut
# nbhd_width: The "width" or "size" or "standard deviation" of the neighbourhood. The neighbourhood weight is a normal distribution with mean of the maternal location & sd nbhd_width
# env_length: The length of the "good" habitat location. This should be varied (i.e. the neighbourhood size to environment size varies) as this may have consequences for dipsersal evo
# t_max: the max number of generations to run the simulation for
# env_change_speed: The shift in the location of the good habitat every time step
# init_loc_mean: Mean of the initial location distribution which initial population locations are drawn from 
# k: From Phillips 2015 - controls the drop in fitness as zw deviates from zero (optimal)
# disp_a_allele: The initial allelic value for the disp_a trait
# disp_b_allele: The initial allelic value for the disp_b trait
# env_allele: The initial allelic value for the environmental trait

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
# t: the generation
# nrows: the number of entries budgetted for the new generation (there is no condition for if the new gen exceeds this number -- should further modify)
# all other input variables as above

make_popn_dataframe <- function(t, nrows, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci){
	
	#empty dataframe with nrow rows and meta_cols+total_genome_length columns
	current_population <- as.data.frame( matrix(data = 0, nrow = nrows, ncol = (meta_cols + total_genome_length)) )

	# give columns the right names
	colnames(current_population) <- meta_col_names
	
	# define some necessary objects - these are used to figure out which columns to create and what information they hold
 	disp_a_locus_1 <- meta_cols+1
	disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
	disp_b_locus_1 <- disp_a_locus_last + 1
	disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
	env_locus_1 <- disp_b_locus_last + 1
	env_locus_last <- env_locus_1 + env_loci*ploidy - 1
	neut_locus_1 <- env_locus_last + 1
	neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

	# Names the columns properly - accounts for changes in number of loci or ploidy of loci controlling traits, as well as differences in number of metadata columns. 
	for (i in (meta_cols+1):ncol(current_population)){
		
		# if the number of metadata columns is even (this matters for deciding which columns correspond to diploid loci pairs, if meta_cols is even then 1st copy is also an even col)
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
			} else {
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
		} else {
		
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
			} else {
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
	
	return(current_population) # empty data frame with the right rows, columns and column labels
}

# This uses the make_popn_dataframe function to create a data frame and then fill it with N individuals
# Initial location of individuals is determined by random draws from norm(init_location, nbhd_width). 
# N: number of inidividuals
# all other input variables as above

make_pop <- function(t, N, init_location, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci){
	
	# make an initial empty data frame (this has one row by default, which is overwritten in the for loop below)
	curr_pop <- make_popn_dataframe(t, 1, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	
	disp_a_locus_1 <- meta_cols+1
	disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
	disp_b_locus_1 <- disp_a_locus_last + 1
	disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
	env_locus_1 <- disp_b_locus_last + 1
	env_locus_last <- env_locus_1 + env_loci*ploidy - 1
	neut_locus_1 <- env_locus_last + 1
	neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

	# fills in the columns with the initial values of the three traits (disp_a, disp_b, and env). Could make random draws from a distribution to seed init pop with genetic variation. 
	
	disp_a_genome <- rep(disp_a_allele/(ploidy*disp_a_loci), ploidy*disp_a_loci)
	disp_b_genome <- rep(disp_b_allele/(ploidy*disp_b_loci), ploidy*disp_b_loci)
	env_genome <- rep(env_allele/(ploidy*env_loci), ploidy*env_loci)
	neut_genome <- rep(0, ploidy*neut_loci)
	
	curr_pop[1:N,] <- 0
	# choose random initial location for every individual in the population, and then create the inidividual. The first zero denotes ???? COME BACK TO THIS
	curr_pop$generation <- t
	curr_pop$individual_ID <- c(1:N)
	curr_pop$location <- rnorm(N, mean = init_location, sd = nbhd_width)
	curr_pop[,disp_a_locus_1:disp_a_locus_last] <- disp_a_genome
	curr_pop[,disp_b_locus_1:disp_b_locus_last] <- disp_b_genome
	curr_pop[,env_locus_1:env_locus_last] <- env_genome
	curr_pop[,neut_locus_1:neut_locus_last] <- neut_genome
	
	return(curr_pop)
}

# -------------------------------- (4) REPRODUCTIVE FUNCTIONS -------------------------------- 

# Calculates the fitness phenotype for an individual from the allelic values at all loci. 
zw <- function(individual, start_locus, end_locus){
	zfit <- sum(individual[start_locus:end_locus])
	zfit <- as.numeric(zfit)
	return(zfit)
}

# Calculates the expected fecundity of an individual (from Phillips 2015)
R <- function(Rmax, zmax, k, zw, nix){
	Ri <- (nix/(nix+0.5))*Rmax*exp(-k*((zw-zmax)^2))
	return(Ri)
}

# Calculates the expected number of offspring (from Phillips 2015 and Beverton & Holt 1958)
expoffspring <- function(zwi, nix, nstar, Rmax, zmax, k){
	Ri <- R(Rmax, zmax, k, zwi, nix)
	alpha <- (Rmax - 1)/nstar
	EO <- Ri/(1+alpha*nix)
	return(EO)
}

# Creates the actual number of offspring for a mother from its expected value (Poisson distributed around a mean of the expected)
reproduce <- function(mom, nstar, Rmax, zmax, k, current_population){
	# nstar is the population size that would be achieved if all individuals had Rmax fecundity. 
	# nix is the current local population density around mom
	
	nix <- localdensity(mom, current_population)

	zwi <- as.numeric(zw(mom, env_locus_1, env_locus_last))

	exp_offspring <- expoffspring(zwi, nix, nstar, Rmax, zmax, k)

	num_offspring <- rpois(1, exp_offspring)
	
	return(num_offspring)	
}

matefinder1D <- function(nbabies, mom, current_population, nbhd_width){

	momID <- mom$individual_ID 
	mom_loc <- mom$location

	myZeroVec <- rep(0, (length(current_population$location)-1)) # SMF comment: the "-1" in the previous line is for taking out the momID
	
	dadIDs <- as.data.frame(cbind(current_population$individual_ID[-momID], current_population$location[-momID], myZeroVec, myZeroVec))
	
	colnames(dadIDs) <- c('individual_ID','location','probability','offspring_counts')
	
	# fills in the probability of mating with each dad given the distance away
	dadIDs$probability <- dnorm(dadIDs$location, mom_loc, nbhd_width) # SMF comment: this is the vectorized step
	
	# list of probs for each dad. Parameterizes a multinomial distribution. Draw from it to get number children each fathers with the given mom
	dadIDs$offspring_counts <- rmultinom(1,size = nbabies, prob=dadIDs$probability)
	
	return(dadIDs)	
}

# this converts the list of dads and their offspring to a format easily fed into the reproduce() function.
convert_dads_list <- function(dadIDs){

	# find the rows in dadIDs where the number of offspring (i.e. column 4) is >0
	dad_indices <- which(dadIDs$offspring_counts>0, arr.ind=TRUE)[,1] 
	
	# find the individual IDs at those rows
	IDs <- as.vector(dadIDs$individual_ID)[dad_indices]
	
	# find the number of offspring (i.e. the value in column 4) for each dad
	offspring <- as.vector(dadIDs$offspring_counts)[dad_indices]
	
	tot_n_off <- sum(offspring)
	
	dads_by_offnum <- vector(length=tot_n_off)
	
	if (length(IDs) > 0){
		placecounter <- 1
		for (i in 1:length(IDs)){
			places <- placecounter:(placecounter+offspring[i]-1)
			dads_by_offnum[places] <- IDs[i]
			placecounter <- placecounter + offspring[i]
		}	
	}

	return(dads_by_offnum)
}

# A function to produce one child to a given mom. This is a helper function for reproduce()
# Assumes that the first two columns of any indidividual have their generation number and their individual number, and col 3 had the locations. So locus 1 is in the fourth entry of the vector that defines every individual, and so on.
# This function cannot yet handle adding new metadata columns!! NEED TO CHANGE THIS STILL


make_offspring <- function(mom, dad, generation, indiv_num, p_mut){ 

	# give the baby it's generation number and individual number
	baby_vec_length = 3 + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1:3] <- c(generation, indiv_num, disperse1D(mom[,3], zda(mom, disp_a_locus_1, disp_a_locus_last), zdb(mom, disp_b_locus_1, disp_b_locus_last)))
	
	# create the baby's genome
	locus_vec <- seq(from = 4, to = total_genome_length+3)
	
	for (i in locus_vec){
		if (i%%2 != 0){ # if this is an the first chromosome of a pair, inherit from mom at random
			momallele = sample(c(i-1,i),1) # choose either the allele at locus i_1 or at locus i_2
			baby[,i] = mom[,momallele]
		}	else { # if this is an odd chromosome, inherit from dad at random
				dadallele = sample(c(i,i+1),1) # choose either the allele at locus i_1 or at locus i_2 from dad's genome
				baby[,i] = dad[,dadallele]	
			}
	}
	
	# now mutate if needed - assumes only one mutation occurs in any individual. This is an okay assumption when the probability of mutating is low and the pop size is not too large.
	if (runif(1) < p_mut){
		focal_locus <- sample(c(4:total_genome_length+3),1) # sample a random locus
		mut_allele <- rnorm(1, baby[,focal_locus], sigma_mut)
		baby[,focal_locus] <- mut_allele
	}
	return(baby)
}



make_offspring2 <- function(mom, dad, generation, indiv_num, p_mut){ 

	# give the baby it's generation number and individual number
	baby_vec_length = 3 + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1:3] <- c(generation, indiv_num, disperse1D(mom[,3], zda(mom, disp_a_locus_1, disp_a_locus_last), zdb(mom, disp_b_locus_1, disp_b_locus_last)))
	
	# create the baby's genome
	locus_vec <- seq(from = 4, to = total_genome_length+3)
	
	for (i in locus_vec){
		if (i%%2 != 0){ # if this is an the first chromosome of a pair, inherit from mom at random
			momallele = sample(c(i-1,i),1) # choose either the allele at locus i_1 or at locus i_2
			baby[,i] = mom[,momallele]
		}	else { # if this is an odd chromosome, inherit from dad at random
				dadallele = sample(c(i,i+1),1) # choose either the allele at locus i_1 or at locus i_2 from dad's genome
				baby[,i] = dad[,dadallele]	
			}
	}
	
	# now mutate if needed - assumes only one mutation occurs in any individual. This is an okay assumption when the probability of mutating is low and the pop size is not too large.
	if (runif(1) < p_mut){
		focal_locus <- sample(c(4:total_genome_length+3),1) # sample a random locus
		mut_allele <- rnorm(1, baby[focal_locus], sigma_mut)
		baby[focal_locus] <- mut_allele
	}
	return(baby)
}

# every individual contributes to the local density based on their distance only -- so here we assume that the fitness phenotype does NOT effect the density that individual contributes. In other words, the every individual's effective density contribution is the same, regardless of its phenotype. 
# This also includes the focal individual in the calculations (this shouldn't matter as long as it's consistent with nstar)

localdensity <- function(mom, current_population){
	
	normdistweights <- sqrt(2*pi*(nbhd_width^2))*dnorm(current_population[,3],mom[,3],nbhd_width)
	local_density <- sum(normdistweights) 

	return(local_density)
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

# A function to choose how far away a given seed disperses.

# -------------------------------- (6) THE ENVIRONMENT -------------------------------- 

# this function provides a map of what the Rmax is at any given point in space asa function of time
# we use the value of t as the middle of the good environment. So at t=0, the good environment extends from [-env_length/2, env_length/2].

environment <- function(dix, Rmax_green, Rmax_red, t, env_length, env_change_speed){ # CLV: 051517 need to update functions that use this
	
	low_lim <- env_change_speed*t-(env_length/2)
	upp_lim <- env_change_speed*t+(env_length/2)
	
	if (low_lim < dix && upp_lim > dix){
		return(c(Rmax_green, zmax_green))
	} else {
		return(c(Rmax_red, zmax_red))
	}
}

# ------------------------------- (7) OTHER FUNCTIONS -----------------------------------

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


