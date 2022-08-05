library("dplyr")
library("tidyverse")
library("ape")
library("parallelsugar")

# function for making a random ancestral sequence
# i ended up just using all "A" for the ancestral sequences
sim_sequence <- function(seq_length, perc_t = 0){
  
  seq_out <- sample(c("A","T"), seq_length, replace = TRUE, 
                    prob = c(1 - perc_t, perc_t))
  
}

# perform mutation on a single base (always A->T, no back mutation)
mutate_base <- function(base, mut_rate){
  
  rand_num <- runif(1) 
  
  if(rand_num <= mut_rate){
    
    base <- "T"
    
  }
  
  return(base)
  
}

# perform mutations across a sequence
mutate_sequence <- function(anc_sequence, mut_rate){
  
  unlist(lapply(anc_sequence, mutate_base, mut_rate = mut_rate))
  
}

# make many sequences from the ancestral sequence via mutation
make_sequences <- function(anc_sequence, n_sequences, mut_rate){
  
  replicate(n_sequences, mutate_sequence(anc_sequence, mut_rate = mut_rate)) %>% 
    data.frame
  
}

# add a fixed proportion of missing genotypes to a site
inject_missing_genotypes <- function(x, prop_missing){
  
  
  n_sites <- floor(length(x)*prop_missing)
  x[sample(1:length(x), n_sites)] <- "N"
  x
  
}

# add missing genotypes across all sites
add_missing_genotypes <- function(seq_df, prop_missing){
  
  
  lapply(seq_df, inject_missing_genotypes, prop_missing = prop_missing) %>%
    data.frame()
  
}

# add missing sites ("N" in all samples) across all sites 
add_missing_sites <- function(seq_df, prop_missing){
  
  n_sites <- floor(nrow(seq_df)*prop_missing)
  sites_to_remove <- sample(1:nrow(seq_df), n_sites)
  empty_row <- seq_df[1,]
  empty_row[] <- "N"
  
  if(n_sites > 0){
    seq_df[sites_to_remove,] <- empty_row
    seq_df
  } else{
    seq_df
  }
  
  
}

# master sequence making funciton
# creates data, labels it, and writes to fasta
make_data_fasta <- function(rep_id, seq_length, mut_rate, n_samples, 
                      prop_missing_genotypes, prop_missing_sites){
  
  anc_sequence <- sim_sequence(seq_length = seq_length)
  seq_df <- make_sequences(anc_sequence, n_sequences = n_samples, mut_rate = mut_rate) %>% 
    add_missing_genotypes(prop_missing = prop_missing_genotypes) %>%
    add_missing_sites(prop_missing = prop_missing_sites)
  
  names(seq_df) <- names(seq_df) %>% 
    paste0("sim_id=", rep_id, ",", 
           "sample=", ., ",",
           "p_missing_genotypes=", prop_missing_genotypes, ",",
           "mut_rate=", mut_rate)
  
  data <- seq_df %>% as.list() 
  
  filename <- paste0("data/fasta/sim=", rep_id, "-", 
                     "p_missing_genotypes=", prop_missing_genotypes, "-",
                     "mut_rate=", mut_rate,
                     ".fasta")
  
  for (i in 1:length(data)){
    
    header <- paste0(">", names(data[i]))
    fasta_concat <- paste0(data[[i]], collapse = "")
    
    write_lines(header, filename, append = TRUE)
    write_lines(fasta_concat, filename, append = TRUE)
    
  }

  
}

# define parameter levels for missing genotypes and sites 
n_missing_genotypes_levels <- seq(0, 0.8, by = 0.1)
mut_rate_levels <- c(0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01)
params <- expand.grid(n_missing_genotypes_levels, mut_rate_levels)

# number of fastas per parameter combination
replicates <- 1000

for (i in 1:nrow(params)){
  
  p_missing_genotypes <- as.numeric(params[i,][1])
  mut_rate <- as.numeric(params[i,][2])
  
  mclapply(1:replicates, make_data_fasta, seq_length = 1000, mut_rate = mut_rate, n_samples = 20, 
                    prop_missing_genotypes = p_missing_genotypes, 
                    prop_missing_sites = 0, mc.cores = 6)
    
}
  




