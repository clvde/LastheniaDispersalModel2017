# Beginning of Dispersal Simulations
# Author: Courtney Van Den Elzen
# Rotation 3 Project: Melbourne and Flaxman Labs
# February 2017 - May 2017

# The goal of this project is to set up a simulation modelling framework which can be used to
# describe the population dynamics and spread of Lasthenia californica (Family Asteraceae), 
# an annual wildflower endemic to california. 

# Assumptions:
	# Discrete generations
	# self-incompatible
	# diploid organisms
	# 5 loci controlling dispersal evolution (there are multiple traits that contribute to the dispersal phenotype)
	# 5 loci controlling environmental evolution (there are multiple traits that contribute to the dispersal phenotype)
	# 5 neutral loci as controls - allow us to show surfing, spatial sorting, etc.
	# no new mutations
	# all seeds from a given mom are fathered by the same dad (this is probably unrealistic but we can change it after)
	# hermaphroditic individuals
	
	
# Structure of the data frame which holds all of the population information
   # Column 1: generation
   # Column 2: individual ID (number from 1 to pop size in the current generation)
   # Column 3-12: dispersal loci (locus 1-1, locus 1-2, locus 2-1, locus 2-2, locus 3-1, locus 3-2, ..., locus 15-1, locus 15-2)
   # Column 13-22: environmental loci
   # Column 23-32: neutral loci
   
   
# Constants
meta_cols <- 3 # the number of metadata columns in the matrix (generation number, individual number, location)
diploid <- 2
disp_loci <- 5
env_loci <- 5
neut_loci <- 5
total_genome_length <- diploid*(disp_loci+env_loci+neut_loci)
max_gens <- 10000
Rmax <- 0 #CHANGE ME

dispersal_kernel_distribution <- #XXXXX CHOOSE THIS BASED ON LITERATURE



# A function to choose dad. This is a helper function for the offspring function

matefinder < - function(current_population, neighbourhood_dist){
	# neighbourhood_dist should be a probability density function that describes the width of a neighbourhood
}

disperse <- function(mom_loc, kernel_dist, kernel_width){
	
	location <- kernel_dist()
}

# A function to produce one child to a given mom. This is a helper function for XXXXX
	
offspring <- function(mom, dad, current_population, generation, indiv_num){ # this whole function assumes that the first two columns of any indidividual have their generation number and their individual number. So locus 1 is in the third entry of the vector that defines every individual, and so on.
	
	# give the baby it's generation number and individual number
	baby_vec_length = metadata + total_genome_length
	baby <- matrix(0,nrow = 1, ncol = baby_vec_length)
	baby[,1] <- generation
	baby[,2] <- indiv_num
	baby[,3] <- disperse()
	
	# create the baby's genome from 
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
	
	
}

# lasthenia californica will produce 1-3 inflorescences total, with 15-30 seeds per inflorescence. 

# A function to describe the environment at every location

# A function to describe the population growth in a given neighbourhood

# A function to choose a mate from the population - weighs probability of mating by distance away

# A function to calculate the dispersal phenotype 

# A function to calculate the fitness phenotype

# A function to calculate dispersal
