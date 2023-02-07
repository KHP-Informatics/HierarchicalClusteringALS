# install packages
BiocManager::install("RNAAgeCalc")
library(RNAAgeCalc)

# define arguments 
expression_data <- args[1] 
tissue_type <- args[2] 
population <- args[3] 
out_dir <- args[4] 

# load in expression data
exprdata <- read.table(expression_data)

# run the program and save the text file
res = predict_age(exprdata = exprdata, tissue = tissue_type, exprtype = "counts", idtype = "ENSEMBL", stype = population, chronage = chronage)
write.table(res, sprintf("%s/RNAAge.txt", out_dir))
