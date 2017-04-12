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
nstar <- 1000
p_mut <- 0.00001 # probability of mutation at every allele
sigma_mut <- 0.001 # standard deviation of the mutation kernel
nbhd_width <- 1 # can set this equal to 1 without loss of generality as the whole environment can be large or small relative to this. 
env_length <- 10 # this should be varied so the gene flow (i.e. the neighbourhood size to environment size varies) as this may have consequences for dipsersal evo
t_max <- 10
init_loc_mean <- 0
nstar <- 10
k <- 1
disp_a_allele <- 1 
disp_b_allele <- 2
env_allele <- 3


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
curr_pop_eg <- make_pop(0, nstar, init_loc_mean, nbhd_width, disp_a_allele, disp_b_allele, env_allele, meta_cols, meta_col_names, ploidy, disp_a_loci, disp_b_loci, env_loci, neut_loci)
mom_eg <- curr_pop_eg[1,]

zwi <- zw(mom_eg, env_locus_1, env_locus_last)
zwi
Rmax <- environment(mom_eg$location, Rmax_good, Rmax_bad, 0, env_length)
Rmax
Ri <- R(Rmax, k, zwi)
Ri
nix <- localdensity(mom_eg, curr_pop_eg)
nix
Eo <- expoffspring(zw, nix, nstar, Rmax_good, k)
nbabies <- reproduce(mom_eg, nstar, Rmax, k, curr_pop_eg)
dadIDs <- matefinder1D(nbabies, mom_eg, curr_pop_eg, nbhd_width)
dadIDs_conv <- convert_dads_list(dadIDs)
zdai <- zda(mom_eg, disp_a_locus_1, disp_a_locus_last)
zdbi <- zdb(mom_eg,disp_b_locus_1, disp_b_locus_last)
baby <- make_offspring(mom_eg, curr_pop_eg[4], 0, 1)
d1D <- disperse1D(mom_eg$location, zdai, zdbi)
d1D
environment(dix, Rmax_good, Rmax_bad, t, env_length)



