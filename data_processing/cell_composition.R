# install packages
install.packages("BRETIGEA")

# load data
expression_data <- args[1]
out_dir <- args[2]
own_marker_expression <- read.table(expression_data)

# run analysis and save data
ct_res = brainCells(own_marker_expression, nMarker = 50, species="human", celltypes = c("ast", "end", "mic", "neu", "oli", "opc"), method="SVD", scale=TRUE)
write.table(ct_res, sprintf("%s/BRETIGEA_celltypeanalysis.txt", out_dir))
