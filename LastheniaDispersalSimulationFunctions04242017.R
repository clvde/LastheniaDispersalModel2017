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
	#empty dataframe
	current_population <- as.data.frame( matrix(data = 0, nrow = nrows, ncol = (meta_cols + total_genome_length)) ) # SMF comment: initially dims were 1 x 1
	# SMF comment: if the data in this object are always going to be zeros, why not just keep it as a matrix? 
	# SMF comment: A matrix has fewer options and properties, so it may not work elsewhere, but just something to think about
	# create the metadata columns
	current_population[,1:meta_cols] <- c(0) # SMF comment: was growing dynamically in columns (and only had one row still)
	# give them the right names
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

	# # create all of the necessary new columns for the loci (after metadata columns) # SMF comment: no longer needed
	# for (i in seq(from=meta_cols, to=(total_genome_length+meta_cols), by=1)){ # SMF comment: grows columns again dynamically; still 1 row
	# 	current_population[,i] <- 0
	# }

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
			} else { # SMF comment: changed to else
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
	
	# add all rows specified above. All entries are zero until changed.
	# current_population[c(1:nrows),] <- 0 # SMF comment: no longer necessary; dynamically grew rows
	return(current_population)
}


# This uses the make_popn_dataframe function to create a data frame and then fill it with N individuals
# Initial location of individuals is determined by random draws from norm(init_location, nbhd_width). 
# N: number of inidividuals
# all other input variables as above
make_pop <- function(t, N, init_location, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci){
	
	# make an initial empty data frame (this has one row by default, which is overwritten in the for loop below)
	curr_pop <- make_popn_dataframe(t, 0, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	
	# fills in the columns with the initial values of the three traits (disp_a, disp_b, and env). Could make random draws from a distribution to seed init pop with genetic variation. 
	# SMF comment: eliminated unnecessary coercion steps
	disp_a_genome <- rep(disp_a_allele/(ploidy*disp_a_loci), ploidy*disp_a_loci)
	disp_b_genome <- rep(disp_b_allele/(ploidy*disp_b_loci), ploidy*disp_b_loci)
	env_genome <- rep(env_allele/(ploidy*env_loci), ncol = ploidy*env_loci)
	neut_genome <- rep(0, ploidy*neut_loci)
	
	# choose random initial location for every individual in the population, and then create the inidividual. The first zero denotes ???? COME BACK TO THIS
	for (i in 1:N){
		init_loc_rand <- rnorm(1, mean = init_location, sd = nbhd_width)
		curr_pop[i,] <- c(0, i, init_loc_rand, disp_a_genome, disp_b_genome, env_genome, neut_genome)
	}
	
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
R <- function(Rmax, k, zw){
	Ri <- Rmax*exp(-k*(zw^2))
	return(Ri)
}

# Calculates the expected number of offspring (from Phillips 2015, but sketchy. Need to validate the reasoning behind this.)
expoffspring <- function(zwi, nix, nstar, Rmax, k){
	Ri <- R(Rmax,k,zwi)
	alpha <- (Rmax - 1)/nstar
	EO <- Ri/(1+alpha*nix)
	return(EO)
}

# Choose pollen donors of offspring. This function does not determine the number of offspring, but instead accepts that as a parameter. 
matefinder1D <- function(nbabies, mom, current_population, nbhd_width){
	# construct a probability density function (multinomial) of mating with every dad (by distance) and then sample from that based on mom's density-independent fecundtiy
	
	# Get mom's identity and location
	momID <- mom$individual_ID
	mom_loc <- mom$location
	
	# Then take mom out of the population (self-incompatibility). This makes the row number different from the individual ID number
	# dads <- current_population[-momID,] # SMF comment: this makes a large new data frame
 	
	# make a data frame of all possible dads with their locations. The third column will hold the probability they father childten with the given mother (determined by distance)
  # SMF comment: the following 3 lines of code, now commented out, were growing the dadIDs object dynamically
	# dadIDs <- as.data.frame(dads$individual_ID)
	# dadIDs[,2] <- dads$location
	# dadIDs[,3] <- 0	
	myZeroVec <- rep(0, (length(current_population$location)-1))
	# SMF comment: the "-1" in the previous line is for taking out the momID
	dadIDs <- as.data.frame(cbind(current_population$individual_ID[-momID], current_population$location[-momID], myZeroVec, myZeroVec))
	colnames(dadIDs) <- c('individual_ID','location','probability','offspring_counts')
	
	# fills in the probability of mating with each dad given the distance away
	# SMF comment: commented out the following loop because it could be vectorized
	# for (i in 1:nrow(dadIDs)){
	# 	prob <- dnorm(dadIDs$location[i], mom_loc, nbhd_width)
	# 	dadIDs[i,3] <- prob
	# }
	dadIDs[,3] <- dnorm(dadIDs$location, mom_loc, nbhd_width) # SMF comment: this is the vectorized step
	
	
	# now have a list of probabilities for each dad, which parameterizes a multinomial distribution. Draw from it to determine the number of children each fathers with the given mom
	# SMF comment: the following two lines were combined into one
	# babies_vec <- rmultinom(1,size = nbabies, prob=dadIDs[,3])
	# dadIDs[,4] <- babies_vec # SMF comment: this was previously a dynamic growth step
	dadIDs[,4] <- rmultinom(1,size = nbabies, prob=dadIDs[,3])
	
	return(dadIDs)	
}


# this converts the list of dads and their offspring to a format easily fed into the reproduce() function.
convert_dads_list <- function(dadIDs){

	dad_indices <- which(dadIDs$offspring_counts>0, arr.ind=TRUE)[,1]
	IDs <- as.vector(dadIDs$individual_ID)[dad_indices]
	offspring <- as.vector(dadIDs$offspring_counts)[dad_indices]
	dads_by_offnum <- vector()
	
	for (i in 1:length(IDs)){
		if (length(offspring) > 0){
			tempvec <- matrix(data = IDs[i], nrow = 1, ncol = offspring[i])
			dads_by_offnum <- c(dads_by_offnum, tempvec)
		}
	}
	
	return(dads_by_offnum)
}

# A function to produce one child to a given mom. This is a helper function for reproduce()
# Assumes that the first two columns of any indidividual have their generation number and their individual number, and col 3 had the locations. So locus 1 is in the fourth entry of the vector that defines every individual, and so on.
# This function cannot yet handle adding new metadata columns!! NEED TO CHANGE THIS STILL

make_offspring <- function(mom, dad, generation, indiv_num){ 

	# give the baby it's generation number and individual number
	baby_vec_length = 3 + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1] <- generation
	baby[,2] <- indiv_num
	baby[,3] <- disperse1D(mom[,3], zda(mom, disp_a_locus_1, disp_a_locus_last), zdb(mom, disp_b_locus_1, disp_b_locus_last))
	
	# create the baby's genome
	locus_vec <- seq(from = 4, to = total_genome_length+3, by = 2)
	for (i in locus_vec){
		if (i%%2 != 0){ # if this is an the first chromosome of a pair, inherit from mom at random
			momallele = round(runif(1)) + i # choose either the allele at locus i_1 or at locus i_2
			baby[,i] = mom[,momallele] 
		}	else { # if this is an odd chromosome, inherit from dad at random
				dadallele = round(runif(1))+i	# choose either the allele at locus i_1 or at locus i_2 from dad's genome
				baby[,i+1] = dad[,dadallele]	
			}
	}
	
	# now mutate if needed - assumes only one mutation occurs in any individual. This is an okay assumption when the probability of mutating is low and the pop size is not too large.
	if (runif(1) < 1){
		focal_locus <- sample(c(4:total_genome_length+3),1) # sample a random locus
		mut_allele <- rnorm(1, baby[focal_locus], sigma_mut)
		baby[focal_locus] <- mut_allele
	}
	return(baby)
}

# every individual contributes to the local density based on their distance only -- so here we assume that the fitness phenotype does NOT effect the density that individual contributes. In other words, the every individual's effective density contribution is the same, regardless of its phenotype. 
# This also includes the focal individual in the calculations (this shouldn't matter as long as it's consistent with nstar)
localdensity <- function(mom, current_population){
	
	distances <- abs(current_population[,3]-mom[3])
	local_density <- 0 
	
	for (i in 1:nrow(current_population)){
		focal_indiv <- current_population[i,]
		dix <- abs(mom$location - focal_indiv$location)
		focal_indiv_contribution <- dnorm(dix, 0, nbhd_width)
		local_density <- local_density + focal_indiv_contribution
	}
	
	return(local_density)
}

# Creates the actual number of offspring for a mother from its expected value (Poisson distributed around a mean of the expected)
reproduce <- function(mom, nstar, Rmax, k, current_population){
	# nstar is the population size that would be achieved if all individuals had Rmax fecundity. 
	# nix is the current local population density around mom
	
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
environment <- function(dix, Rmax_good, Rmax_bad, t, env_length, env_change_speed){
	
	low_lim <- env_change_speed*t-(env_length/2)
	upp_lim <- env_change_speed*t+(env_length/2)
	
	if (low_lim < dix && upp_lim > dix){
		return(Rmax_good)
	} else {
		return(Rmax_bad)
	}
}


