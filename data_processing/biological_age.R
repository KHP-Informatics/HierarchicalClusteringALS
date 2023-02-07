# this script works after you have git cloned https://github.com/gemmashireby/CorticalClock

# arguments
args[1] = clock_directory
args[2] = beta_matrix # rownames=cpgs, colnames=IDs
args[3] = phenotype_file
args[4] = phenotype_ids
args[5] = phenotype_ages
args[6] = out_dir

# source the cortical clock script
source(sprintf("%s/PredCorticalAge/CorticalClock.r", clock_directory))

# run the script
setwd(out_dir)
CorticalClock(beta_matrix, phenotype_file, sprintf("%s/PredCorticalAge", clock_directory), phenotype_file$phenotype_ids, phenotype_file$phenotype_ages) 
