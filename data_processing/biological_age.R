# this script works after you have git cloned https://github.com/gemmashireby/CorticalClock

# arguments
clock_directory <- args[1]
beta_matrix <- args[2]  # rownames=cpgs, colnames=IDs
phenotype_file <- args[3] 
phenotype_ids <- args[4]
phenotype_ages <- args[5] 
out_dir <- args[6] 

# source the cortical clock script
source(sprintf("%s/PredCorticalAge/CorticalClock.r", clock_directory))

# run the script
setwd(out_dir)
CorticalClock(beta_matrix, phenotype_file, sprintf("%s/PredCorticalAge", clock_directory), phenotype_file$phenotype_ids, phenotype_file$phenotype_ages) 
