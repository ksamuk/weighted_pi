library("tidyverse")
library("ape")
library("pegas")
library("parallelsugar")

# original PIs function
PIs <- function(data, ...){
  dnas <- as.data.frame(as.character(data))
  freqs <- apply(dnas,1,table,exclude = "N")
  nDiff <- sapply(freqs,function(z){ 
    tabs <- outer(z,z)
    sum(tabs[lower.tri(tabs)])})
  nComp <- sapply(freqs,sum)
  nComp <- (nComp*(nComp-1))/2
  WeightedPi <- mean(nDiff/nComp,na.rm = TRUE)
  PEGASres <- pegas::nuc.div(data, pairwise_deletion = TRUE)
  PIXYres <- sum(nDiff)/sum(nComp)
  return(data.frame(weighted = WeightedPi,
                    pixy = PIXYres, pegas = PEGASres))
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



# list all fasta files to be analyzed
file_names <- list.files("data/fasta", full.names = TRUE)

# function for reading in file, parsing parameters, and computing pis
read_data_compute_pis <- function(file_name, pwd_toggle){
  
  #print(file_name)
  data <- read.FASTA(file_name)
  p_missing_genotypes <- gsub(".*_genotypes=|-mut_rate.*", "", file_name)
  mut_rate <- gsub(".*-mut_rate=|\\.fas.*", "", file_name)
  data.frame(file_name, p_missing_genotypes, mut_rate, PIs_fixed(data, pwd_toggle))
  
}

# perform computations and bind to data frame
results_df <- mclapply(file_names, read_data_compute_pis, pwd_toggle = TRUE, mc.cores = 6) %>%
  bind_rows() 

results_df <- results_df %>%
  type_convert()

# save for recovery elsewhere
saveRDS(results_df, file = "data/simulation_data_results.rds")

results_df <- read_rds("data/simulation_data_results.rds")

# no evidence of change in median value across multiple levels of diversity
# and levels of missing data
figure2 <- results_df %>%
  mutate(exp_pi = signif(1 - (((1-mut_rate)^2)+(mut_rate^2)), 2)) %>%
  mutate(exp_pi = ifelse(exp_pi == 0.0099, 0.01, exp_pi)) %>% #corrects rounding error from about
  gather(key = "method", value = "estimate", -file_name, -p_missing_genotypes, -mut_rate, -exp_pi) %>%
  mutate(method = ifelse(method == "weighted", "unweighted", method)) %>%
  filter(method != "pegas_pi") %>%
  ggplot(aes(x = method, y = estimate))+
  geom_boxplot(outlier.size=0.5)+
  facet_grid(exp_pi~p_missing_genotypes, scales = "free")+
  xlab("Method")+
  ylab("Estimated average pi")+ 
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Expected pi", breaks = NULL, labels = NULL))

ggsave(filename = "figures/Figure2_raw.pdf", figure2, device = "pdf", width = 8, height = 6)
