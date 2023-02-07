# install packages
BiocManager::install("RNAAgeCalc")
library(RNAAgeCalc)

# define arguments 
args[1] = expression_data
args[2] = tissue_type
args[3] = population
args[4] = out_dir

# load in expression data
exprdata <- read.table(expression_data)

# run the program and save the text file
res = predict_age(exprdata = exprdata, tissue = tissue_type, exprtype = "counts", idtype = "ENSEMBL", stype = population, chronage = chronage)
write.table(res, sprintf("%s/RNAAge.txt", out_dir))
