# Beginning of Dispersal Simulations
# Author: Courtney Van Den Elzen
# Rotation 3 Project: Melbourne and Flaxman Labs
# February 2017 - May 2017

# The goal of this project is to set up a simulation modelling framework which can be used to describe the population dynamics and spread of annual, wind-dispersed grassland plants

# Assumptions:
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
	
	
# Structure of the data frame which holds all of the population information
   # Column 1: generation
   # Column 2: individual ID (number from 1 to pop size in the current generation)
   # Column 3: location
   # Column 4-13: dispersal trait 1 (a = mean) loci (col3 = locus 1-1, col4 = locus 1-2, col5 = locus 2-1, col6 = locus 2-2, col7 = locus 3-1, col8 = locus 3-2, col9 = locus 4-1, col10 = locus 4-2, col11 = locus 5-1, col12 = locus 5-2)
   # Column 14-23: dispersal trait 2 (b = shape) loci 
   # Column 24-33: environmental loci
   # Column 34-43: neutral loci

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
# Constants

meta_cols <- 3 # the number of metadata columns in the matrix (generation number, individual number, location)
diploid <- 2
disp_loci <- 5
env_loci <- 5
neut_loci <- 5
total_genome_length <- diploid*(disp_loci+env_loci+neut_loci)
fit_start_locus <- meta_cols+1
fit_end_locus <- fit_start_locus+diploid*env_loci
disp_start_locus <- meta_cols+1
disp_end_locus 
max_gens <- 10000
Rmax <- 150 # (50 seeds per inflorescence, 3 inflorescences MAX)
dispersal_kernel_distribution <- # WALD (inverse gaussian)
sigma <- 1 # effective size of spatial neighbourhood
nstar <- 1000
alpha <- (Rmax - 1)/nstar
b <- 1 # shape parameter of the inverse gaussian distribution
p_recomb <- 0.0001
p_mut <- 0.00001 # probability of mutation at every allele
sigma_mut <- 0.001 # standard deviation of the mutation kernel
nbhd_width <- 1 # can set this equal to 1 without loss of generality as the whole environment can be large or small relative to this. 
env_length <- ####
env_width <- ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculates the expected fecundity of an individual

R <- function(Rmax,k,zw){
	Ri <- Rmax*exp(-k*(zw^2))
	return(Ri)
}


# Calculates the expected number of offspring (as in Phillips 2015)

expoffspring <- function(Rmax,k,zw,Nix){
	Ri <- R(Rmax,k,zw)
	EO <- Ri/(1+alpha*Nix)
	return(EO)
}


# Calculates the first dispersal phenotype (a- mean distance) for an individual.

zda <- function(individual,start_locus,total_loci){
	# this assumes that you have one vector describing an individual. entries 13-22 describe the dispersal alleles
	end_locus <- total_loci+start_locus-1
	zdispa <- sum(individual[start_locus:end_locus])
	return(zdispa)
}


zdb <- function(individual,start_locus,total_loci){
	# this assumes that you have one vector describing an individual. entries 13-22 describe the dispersal alleles
	end_locus <- total_loci+start_locus-1
	zdisp <- sum(individual[start_locus:end_locus])
	return(zdispa)
}


# Calculates the fitness phenotype for an individual

zw <- function(individuals,start_locus,total_loci){
	end_locus <- total_loci+start_locus-1
	zfit <- sum(individual[start_locus:end_locus])
	return(zfit)
}


# A function to choose dad. This is a helper function for the offspring function


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

disperse2D <- function(mom_loc, zda, zdb){
	
	xm <- momloc[1] # assumes moms location
	ym <- momloc[2]
	theta <- 360*runif(1)
	dist <- rinvgauss(1, zda, zdb)
	
	if (theta == 0 || theta == 360) {
		
		xb <- xm
		yb <- ym + dist
		
	} else if (0 < theta < 90) {
		
		xb <- xm - dist*sin(theta)
		yb <- ym + dist*cos(theta)
		
	} else if (theta == 90) {
		
		xb <- xm - dist
		yb <- ym
		
	} else if (90 < theta < 180) {
		
		xb <- xm - dist*sin(180-theta)
		yb <- ym - dist*cos(180-theta)
		
	} else if (theta == 180) {
		
		xb <- xm 
		yb <- ym - dist
				
	} else if (180 < theta < 270) {
		
		xb <- xm + dist*sin(270-theta)
		yb <- ym - dist*cos(270-theta)
		
	} else if (theta == 270) {
		
		xb <- xm + dist
		yb <- ym
		
	} else if (270 < theta < 360) {
		
		xb <- xm + dist*sin(360-theta)
		yb <- ym + dist*cos(360-theta)
		
	}	
	
	baby_loc <- c(xb,yb)
	
	return(baby_loc)
	
}


matefinder1D < - function(mom, current_population, neighbourhood_dist){
	# neighbourhood_dist should be a probability density function that describes the width of a neighbourhood
	
	# There are two ways I could write this function: (1) randomly choose a location using the normal distribution and then find the nearest individual to that point to mate with OR (2) find the probability of 	mating with every individual and then sample from that discrete distribution. I have a feeling the former will be faster? Especially if I could somehow search nearest that point for individuals. I should 	talk to Brett and Sam about this. 
	
	mate_loc <- rnorm(1,mom[3],nbhd_width)
	popsize <- nrow(current_population)
	locations <- current_population[,3]
	mate_dists <- locations - mate_loc
	min_list <- which(mate_dists == min(mate_dists))
	if (length(min_list) > 1){
		rand_choice <- sample(min_list,1)
		mate <- current_population[rand_choice,]
	} else{
		mate <- min_list[1]
	}
	return(mate)	
}


# A function to produce one child to a given mom. This is a helper function for XXXXX
	
offspring <- function(mom, dad, current_population, generation, indiv_num){ # this whole function assumes that the first two columns of any indidividual have their generation number and their individual number. So locus 1 is in the third entry of the vector that defines every individual, and so on.
	
	# give the baby it's generation number and individual number
	baby_vec_length = metadata + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1] <- generation
	baby[,2] <- indiv_num
	baby[,3] <- disperse1D(mom[,3],zda(mom),zdb(mom))
	
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
	
	
	
	
}

# Creates the drawn number of offspring for a given mother and packages them in a matrix

reproduce <- function(mom, dad, nstar, Nix){
	# nstar is the population size that would be achieved if all individuals had Rmax fecundity. 
	# Nix is the current local population size
	
	
	
}
# lasthenia californica will produce 1-3 inflorescences total, with 15-30 seeds per inflorescence. 

# step 1 - find a mate (matefinder)
# step 2 - reproduce (offspring - makes one offspring, reproduce - calculates the number of offspring and then gets offspring function to make them)
# step 3 - 
