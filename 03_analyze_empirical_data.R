library("tidyverse")
library("ape")
library("pegas")
library("snow")
library("parallel")
#devtools::install_github("https://github.com/nathanvan/parallelsugar")
library("parallelsugar")

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


# fixed PIs function based on pixy code
PIs_fixed <- function(data, pwd_toggle){
  
  dnas <- as.data.frame(as.character(data))
  
  # create table of alleles from FASTA file
  allele_count_tab <- apply(dnas,1,table)
  
  if(is.numeric(allele_count_tab)){
    
    allele_count_tab <- data.frame(a = allele_count_tab)
    
  } else{
    
    allele_count_tab <- allele_count_tab %>% 
      bind_rows() %>% 
      mutate_if(is.table, as.numeric)
    
  }
  
  # remove n counts (not used)
  if ("n" %in% names(allele_count_tab)){
    
    allele_count_tab <- allele_count_tab %>% 
      select(-n)
    
  }
  
  if ("N" %in% names(allele_count_tab)){
    
    allele_count_tab <- allele_count_tab %>% 
      select(-N)
    
  }
  
  # convert alleles to major + minor allele counts
  
  # initialize a new dataframe 
  maj_min_tab <- data.frame(maj = rep(NA, nrow(allele_count_tab)),
                            min = rep(NA, nrow(allele_count_tab)))
  
  
  # assign major/minor allele 
  for (i in 1:nrow(allele_count_tab)){
    
    alleles <- allele_count_tab[i,]
    
    # if there are more than two alleles at the size, drop it
    if (sum(alleles != 0, na.rm = TRUE) > 2){
      
      maj_min_tab[i,] <- data.frame(maj = NA, min = NA)
      
    } else{
      
      alleles_sorted <- alleles %>% as.numeric %>% sort(decreasing = TRUE)
      
      if (length(alleles_sorted) == 1){
        
        maj_min_tab[i,] <- data.frame(maj = alleles_sorted, min = 0)
        
      } else{
        
        maj_min_tab[i,] <- data.frame(maj = alleles_sorted[1], min = alleles_sorted[2])
        
        
      }
      
    }
    
  }
  
  # calculate pi from maj/minor alleles
  pi_calc_df <- maj_min_tab %>%
    mutate(n_diffs = maj * min) %>%
    mutate(n_gts = maj + min) %>%
    mutate(n_comps = choose(n_gts,2)) %>%
    mutate(pi_w_site = n_diffs/n_comps)
  
  weighted_pi <- mean(pi_calc_df$pi_w_site, na.rm = TRUE)
  pixy_pi <- sum(pi_calc_df$n_diffs, na.rm = TRUE)/sum(pi_calc_df$n_comps, na.rm = TRUE)
  
  PEGASres <- pegas::nuc.div(data, pairwise_deletion = pwd_toggle)
  
  return(data.frame(weighted = weighted_pi,
                    pixy = pixy_pi, pegas_pi = PEGASres))
  
}
# results are identical


# confirm that results are identical with no missing data
# cat data
cat_data <- read.FASTA("data/cat.fasta")
PIs_fixed(cat_data)


# barb data
barb_data <- read.FASTA("data/barb.fasta")
PIs_fixed(barb_data)


# WITH MISSING DATA:

# function with maciek's method of adding missing data
# altered to include pi method from pixy code
ftPIs <- function(index, dane){
  
  fasta <- t(data.frame(as.character(dane)))
  xN <- fasta
  
  Ndistr <- NA
  Ndistr1 <- c()
  Ndistr2 <- c()
  
  # this doesn't seem to affect the amount of missing data
  propN <- sample(seq(0.5,0.9,by = 0.1),1)
  
  while(is.na(mean(Ndistr))){
    while(length(c(Ndistr1,Ndistr2))<ncol(xN)){
      
      Ndistr1 <- round(rgamma(round(ncol(xN)*propN), sample(seq(0,10,by = 0.1),1),1))
      Ndistr2 <- round(runif(round(ncol(xN)*(1-propN)),min = 0, max=round(nrow(xN)*0.80)))
      
    }
    Ndistr <- c(Ndistr1,Ndistr2)[sample(1:length(c(Ndistr1,Ndistr2)),ncol(xN))]
  }
  for(j in 1:ncol(xN)){
    if(Ndistr[j]<ncol(xN)){
      xN[sample(1:nrow(xN),Ndistr[j]),j] <- "N"} else xN[j,] <- "N" 
  }
  
  freqs <- apply(xN,2,table,exclude = "N")
  emptySites <- sapply(unlist(freqs),dim)
  
  if (!is.list(emptySites)) xN <- xN[,emptySites!=0]
  
  # FOR MISSING DATA CASE
  allele_count_tab <- apply(xN,2,table) %>% 
    bind_rows() %>% 
    mutate_if(is.table, as.numeric) %>%
    select(-N)
  
  # convert alleles to major + minor allele counts
  
  # initialize a new dataframe 
  maj_min_tab <- data.frame(maj = rep(NA, nrow(allele_count_tab)),
                            min = rep(NA, nrow(allele_count_tab)))
  
  
  # assign major/minor allele 
  for (i in 1:nrow(allele_count_tab)){
    alleles <- allele_count_tab[i,]
    # if there are more than two alleles at the site, drop it
    if (sum(alleles != 0, na.rm = TRUE) > 2){
      maj_min_tab[i,] <- data.frame(maj = NA, min = NA)
    } else{
      alleles_sorted <- alleles %>% as.numeric %>% sort(decreasing = TRUE)
      if (length(alleles_sorted) == 1){
        maj_min_tab[i,] <- data.frame(maj = alleles_sorted, min = 0)
      } else{
        maj_min_tab[i,] <- data.frame(maj = alleles_sorted[1], min = alleles_sorted[2])
      }
    }
  }
  
  # calculate pi from maj/minor alleles
  pi_calc_df <- maj_min_tab %>%
    mutate(n_diffs = maj * min) %>%
    mutate(n_gts = maj + min) %>%
    mutate(n_comps = choose(n_gts,2)) %>%
    mutate(pi_w_site = n_diffs/n_comps) 
  
  
  
  weighted_pi <- mean(pi_calc_df$pi_w_site, na.rm = TRUE)
  pixy_pi <- sum(pi_calc_df$n_diffs, na.rm = TRUE)/sum(pi_calc_df$n_comps, na.rm = TRUE)
  
  PEGASres <- pegas::nuc.div(dane)
  
  return(data.frame(Ns = sum(xN=="N")/length(xN), weighted = weighted_pi,
                    pixy = pixy_pi, pegas_pi = PEGASres))
  
}


barbs <- read.FASTA("data/barb.fasta",type = "DNA")

# original missing distribution
barbusSims <- parallelsugar::mclapply(1:2000, ftPIs, dane = barbs, 
                                      mc.cores = 6)
# save to avoid rerunning
saveRDS(barbusSims, file = "data/barbus_sims_orig_out.rds")

# save to avoid rerunning
barbusSims <- read_rds("data/barbus_sims_orig_out.rds")

# this shows the pattern from the konopinski paper
bind_rows(barbusSims) %>%
  select(-pegas_pi) %>%
  gather(key = "method", value = "estimate", -Ns) %>%
  ggplot(aes(x = estimate, color = method))+
  geom_density()

# rough recreation of the boxplot
# under this specific missing data model, 
# the unweighted estimator appears to do better 
bind_rows(barbusSims) %>%
  select(-pegas_pi) %>%
  mutate(Ns_cut = cut_interval(Ns, 4)) %>%
  select(-Ns) %>%
  gather(key = "method", value = "estimate", -Ns_cut) %>%
  ggplot(aes(x = method, y= estimate))+
  geom_boxplot()+
  facet_wrap(~Ns_cut)

# uniform missing data distribution

ftPIs_uniform <- function(index, dane,  prop_missing_levels){
  
  results_list <- list()
  
  for (j in 1:length(prop_missing_levels)){
    
    fasta <- data.frame(as.character(dane))
    
    prop_missing <- prop_missing_levels[j]
    
    fasta <- fasta %>%
      add_missing_genotypes(prop_missing)
    
    
    # FOR MISSING DATA CASE
    if(prop_missing > 0 ){
      
      allele_count_tab <- apply(fasta,1,table) %>% 
        bind_rows() %>% 
        mutate_if(is.table, as.numeric) %>%
        select(-N)
      
    } else{
      
      allele_count_tab <- apply(fasta,1,table) %>% 
        bind_rows() %>% 
        mutate_if(is.table, as.numeric) 
      
    }

    
    # convert alleles to major + minor allele counts
    
    # initialize a new dataframe 
    maj_min_tab <- data.frame(maj = rep(NA, nrow(allele_count_tab)),
                              min = rep(NA, nrow(allele_count_tab)))
    
    
    # assign major/minor allele 
    for (i in 1:nrow(allele_count_tab)){
      alleles <- allele_count_tab[i,]
      # if there are more than two alleles at the site, drop it
      if (sum(alleles != 0, na.rm = TRUE) > 2){
        maj_min_tab[i,] <- data.frame(maj = NA, min = NA)
      } else{
        alleles_sorted <- alleles %>% as.numeric %>% sort(decreasing = TRUE)
        if (length(alleles_sorted) == 1){
          maj_min_tab[i,] <- data.frame(maj = alleles_sorted, min = 0)
        } else{
          maj_min_tab[i,] <- data.frame(maj = alleles_sorted[1], min = alleles_sorted[2])
        }
      }
    }
    
    # calculate pi from maj/minor alleles
    pi_calc_df <- maj_min_tab %>%
      mutate(n_diffs = maj * min) %>%
      mutate(n_gts = maj + min) %>%
      mutate(n_comps = choose(n_gts,2)) %>%
      mutate(pi_w_site = n_diffs/n_comps) 
    
    weighted_pi <- mean(pi_calc_df$pi_w_site, na.rm = TRUE)
    pixy_pi <- sum(pi_calc_df$n_diffs, na.rm = TRUE)/sum(pi_calc_df$n_comps, na.rm = TRUE)
    
    PEGASres <- pegas::nuc.div(dane)
    
    results_list[[j]] <- data.frame(index, prop_missing = prop_missing, 
                                    weighted = weighted_pi, pixy = pixy_pi, 
                                    pegas_pi = PEGASres)
    
  }
  
  return(results_list)
  
}


barbusSims1 <- parallelsugar::mclapply(1:1000, ftPIs_uniform, dane = barbs, 
                                       prop_missing_levels = seq(from = 0.1, to = 0.5, by = 0.1),
                                       mc.cores = 6)

saveRDS(barbusSims1, file = "data/barbus_sims_uniform_out.rds")

barbusSims1 <- read_rds("data/barbus_sims_uniform_out.rds")

figure1 <- bind_rows(barbusSims1) %>%
  select(-pegas_pi, -index) %>%
  gather(key = "method", value = "estimate", -prop_missing) %>%
  mutate(method = ifelse(method == "weighted", "unweighted", method)) %>%
  ggplot(aes(x = method, y= estimate))+
  geom_boxplot()+
  facet_grid(~prop_missing)

ggsave(filename = "figures/Figure1_raw.pdf", figure1, device = "pdf", width = 5, height = 3)


