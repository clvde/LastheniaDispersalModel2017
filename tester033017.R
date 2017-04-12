meta_cols <- 3 # the number of metadata columns in the matrix (generation number, individual number, location)
diploid <- 2
disp_a_loci <- 5
disp_b_loci <- 5
env_loci <- 5
neut_loci <- 5
total_genome_length <- diploid*(disp_a_loci+disp_b_loci+env_loci+neut_loci)
max_gens <- 10000
Rmax <- 150 # (50 seeds per inflorescence, 3 inflorescences MAX)
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

current_population <- data.frame('generation' = 1, 'individual_ID' = 1:11, 'location' = seq(from = 0, to = 1, by = .1),'rand' = seq(from = 0, to = 1, by = .1))

meta_cols=ncol(current_population)
disp_a_locus_1 <- meta_cols+1
disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*diploid - 1
disp_b_locus_1 <- disp_a_locus_last + 1
disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*diploid - 1
env_locus_1 <- disp_b_locus_last + 1
env_locus_last <- env_locus_1 + env_loci*diploid - 1
neut_locus_1 <- env_locus_last + 1
neut_locus_last <- neut_locus_1 + neut_loci*diploid - 1

# creates all of the necessary new columns for the loci
for (i in seq(from=meta_cols, to=total_genome_length+meta_cols,by=1)){
	current_population[,i] <- 0
}

# Names the columns properly

for (i in meta_cols+1:ncol(current_population)){
	# j = i-meta_cols
	# k = round((j/2)+0.0000001) # round function in r round 0.5 to the even digit (i.e. 0. This is so weird.)
	
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

