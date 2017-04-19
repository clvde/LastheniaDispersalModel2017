#Packages
library(statmod)

# Constants
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
nstar <- 100
p_mut <- 0.00001 # probability of mutation at every allele
sigma_mut <- 0.001 # standard deviation of the mutation kernel
nbhd_width <- 1 # can set this equal to 1 without loss of generality as the whole environment can be large or small relative to this. 
env_length <- 10 # this should be varied so the gene flow (i.e. the neighbourhood size to environment size varies) as this may have consequences for dipsersal evo
t_max <- 10
init_loc_mean <- 0
k <- 1
disp_a_allele <- 1 
disp_b_allele <- 2
env_allele <- 0

# Derived Constants
disp_a_locus_1 <- meta_cols+1
disp_a_locus_last <- disp_a_locus_1 + disp_a_loci*ploidy - 1
disp_b_locus_1 <- disp_a_locus_last + 1
disp_b_locus_last <- disp_b_locus_1 + disp_b_loci*ploidy - 1
env_locus_1 <- disp_b_locus_last + 1
env_locus_last <- env_locus_1 + env_loci*ploidy - 1
neut_locus_1 <- env_locus_last + 1
neut_locus_last <- neut_locus_1 + neut_loci*ploidy - 1

# derived inputs/function tests

df_shell <- make_popn_dataframe(t, 3, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
df_shell
curr_pop_eg <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)

mom_eg <- curr_pop_eg[1,]
curr_pop_eg <- curr_pop_eg[-1,]
# mom_eg

# curr_pop_eg

# zwi <- zw(mom_eg, env_locus_1, env_locus_last)
# zwi
# Rmax <- environment(mom_eg$location, Rmax_good, Rmax_bad, 0, env_length)
# Rmax
# Ri <- R(Rmax, k, zwi)
# Ri
# nix <- localdensity(mom_eg, curr_pop_eg)
# nix
# Eo <- expoffspring(zwi, nix, nstar, Rmax_good, k)
# Eo
# nbabies <- reproduce(mom_eg, nstar, Rmax, k, curr_pop_eg)
# nbabies
# dadIDs <- matefinder1D(nbabies, mom_eg, curr_pop_eg, nbhd_width)
# dadIDs
# dadIDs_conv <- convert_dads_list(dadIDs)
# dadIDs_conv
# zdai <- zda(mom_eg, disp_a_locus_1, disp_a_locus_last)
# zdai
# zdbi <- zdb(mom_eg,disp_b_locus_1, disp_b_locus_last)
# zdbi
# baby <- make_offspring(mom_eg, curr_pop_eg[4], 0, 1)
# baby
# d1D <- disperse1D(mom_eg$location, zdai, zdbi)
# d1D
# environment(mom_eg$location, Rmax_good, Rmax_bad, 0, env_length)

# Now Simulate!
current_population <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)

write_name <- paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Practice Files/lascali_sim_dispevo_only_gen_",0,".csv", sep="")
write.csv(current_population, write_name)

for (t in 1:t_max){
	# (1) Reproduction
	# (2) Parental Death
	# (3) Dispersal (but this is density independent)
	# (3) F1 Reproduction
	
	# (1) & (2) - offspring dispersal is built into the make_offspring function. 
	print(t)
	next_generation <- make_popn_dataframe(t, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
	next_gen_ID_tracker <- 1
	
	for (i in 1:nrow(current_population)){
		mom <- current_population[i,]
		Rmax <- environment(mom$location, Rmax_good, Rmax_bad, t, env_length)
		mates <- current_population[-i,]
		n_offspring <- reproduce(mom, nstar, Rmax, k, mates)
		dads_list <- matefinder1D(n_offspring, mom, mates, nbhd_width)
		dads_list_reformat <- convert_dads_list(dads_list)
		
		if (n_offspring > 0) {
			for (n in 1:n_offspring){
				dad <- dads_list_reformat[n]
				offspring <- make_offspring(mom, dad, t, next_gen_ID_tracker)
				next_generation[next_gen_ID_tracker,] <- offspring
				next_gen_ID_tracker <- next_gen_ID_tracker + 1
			}
		}
	}
	
	# write the parental generation to file before erasing them (annuals)
	write_name <- paste("/Users/Courtney/Documents/Rotation 3 - Melbourne & Flaxman Labs/Simulation Practice Files/lascali_sim_dispevo_only_gen_",t,".csv", sep="")
	current_population <- next_generation
	write.csv(current_population, write_name)
	
	if (nrow(current_population) == 0){
		print("extinction!")
		break
	}
}

