## 1. Install and load libraries
BiocManager::install("biomaRt")
install.packages("Hmisc")
BiocManager::install("DESeq2")
install.packages("NMF")
BiocManager::install("BiocParallel")
install.packages("knitr")
library("biomaRt")
library("Hmisc")
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(20))
library("knitr")
library("NMF")

## 2. Load in the data (cases only) and perform VST normalisation using DESeq2 (no design)
countdata <- read.table("RNAseq_counts_cases.txt", row.names=1, header=TRUE)
head(countdata)
coldata <- read.table("RNAseq_samples_cases_info.txt", row.names=1, header=TRUE)
head(coldata)
dds_genes <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ 1)
dds_genes <- estimateSizeFactors(dds_genes)
idx <- rowSums(counts(dds_genes, normalized=TRUE) >=5 ) >= 10
dds_genes <- dds_genes[idx,]
dds_genes
vsd <- vst(dds_genes, blind=TRUE) #this means it calculates within group variability
head(assay(vsd),3)
vsd.out <- assay(vsd)
write.table(vsd.out,file="RNAseq_cases_VSD.txt",sep="\t")

## 3. Remove sex chromosomes 
ensembl = useEnsembl(biomart="ensembl",GRCh=38, dataset="hsapiens_gene_ensembl")
head(listFilters(ensembl))
head(listAttributes(ensembl))
positions <- read.table("hg38_sex_chrom_positions.txt")
ensembl_sex <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), filters = c("chromosome_name", "start", "end"), values = list(positions[,1], positions[,2], positions[,3]), mart = ensembl)
ensembl_sex_name <- list(ensembl_sex$ensembl_gene_id)

vsd_nosex_list <- setDT(vsd.out)
vsd_nosex <- vsd.out[rownames(vsd.out) %nin% unlist(ensembl_sex_name)]
dim(vsd_nosex)
write.table(vsd_nosex, "RNAseq_cases_VSD_nosex.txt")

## 4. Extract the top 5000 variably expressed genes
row_MADs <- apply(log(vsd_nosex + 0.1), 1, mad)
row_MADs_df <- data.frame(row_MADs)
n_5000_genes <- head(row_MADs_df[order(row_MADs_df$row_MADs, decreasing = TRUE), ,drop=FALSE], 5000)
n_5000_genes_list <- list(rownames(n_5000_genes))

vsd_nosex_top5000 <- vsd_nosex[rownames(vsd_nosex) %in% unlist(n_5000_genes_list), ]
dim(vsd_nosex_top5000)
write.table(vsd_nosex_top5000, "RNAseq_cases_VSD_nosex_top5000.txt")

## 5. Run NMF clustering (several functions taken from the SAKE package https://github.com/naikai/sake)

## 5.1 Select K setting with the highest cophenetic correlation coefficient
estim.r <- nmf(vsd_nosex_top5000, 2:10, nrun=100, "nsNMF", .opt=paste0("vp", 20), seed=123211, maxIter=1000)
res <- estim.r
nmf_rank <- estim.r$measures[, 1]
### plot for cophenetic correlation
ylabel <- colnames(estim.r$measures)[18]
plot(nmf_rank, estim.r$measures[,13], type="o", xlab="Rank", ylab=ylabel)
### plot for dispersion
ylabel <- colnames(estim.r$measures)[19]
plot(nmf_rank, estim.r$measures[,13], type="o", xlab="Rank", ylab=ylabel)
### plot for consensus
ylabel <- colnames(estim.r$measures)[21]
plot(nmf_rank, estim.r$measures[,13], type="o", xlab="Rank", ylab=ylabel)

## 5.2 Run NMF with selected k setting (was 3 in our case)
res2 <- nmf(vsd_nosex_top5000, 3, nrun=100, "nsNMF", .opt=paste0("vp", 20), seed=123211, maxIter=1000)

## 5.3 Generate different summary plots
### summary of measures used to evaluate which K setting is best
summary.class <- as.data.frame(NMF::summary(res2))
write.csv(summary.class, "estimateK1000iterations.txt")

### extract most informative features (genes) which define each of the clusters, and the total gene assignment
nmf_extract_feature <- function(res2, rawdata=NULL, manual.num=0, method="default", math="mad", FScutoff=0.9){
   feature.score <- featureScore(res2)
   predict.feature <- predict(res2, what="features", prob=T)

   data.feature <- data.frame(Gene=names(feature.score),
                              featureScore=feature.score,
                              Group=predict.feature$predict,
                              prob=predict.feature$prob,
                              stringsAsFactors=FALSE)
   if(method=="total"){
      print("return all featureScores")
   }else{
      if(method=="default"){
         # extracted features for each group
         if (manual.num==0){
            extract.feature <- extractFeatures(res2)
         }else if (manual.num>0 && manual.num<=length(featureNames(res2))){
            extract.feature <- extractFeatures(res2, manual.num)
         }else{
               stop("wrong number of (manual num) features ")
         }

         data.feature <- extract.feature %>%
                     lapply(., function(x) data.feature[x, ]) %>%
                     rbindlist %>%
                     as.data.frame.matrix
      }else if(method=="rank"){
         if(is.null(rawdata)){
            stop("error: need to provide original expression data if method is 'rank' ")
         }
         # data.feature <- cbind(data.feature, math=apply(log2(rawdata+1), 1, math)) %>%
         data.feature <- cbind(data.feature, math=apply(rawdata, 1, math)) %>%
                              filter(featureScore>=FScutoff) %>%
                              arrange(Group, dplyr::desc(math), dplyr::desc(prob))
         if(manual.num>0){
            data.feature <- group_by(data.feature, Group) %>%
                              top_n(manual.num)
                              # top_n(manual.num, math)
                              # filter(min_rank(desc(math))<=manual.num)
         }
      }
   }
   return(data.feature)
}
                            
write.table(nmf_extract_feature(res2, method="default", math="mad"), "genes3clusters_mostinformativegenespercluster.txt")
write.table(nmf_extract_feature(res2, method="total", math="mad"), "genes3clusters_totalgenespercluster.txt")

### extract groups i.e. samples within clusters
nmf_extract_group <- function(res, type="consensus", matchConseOrder=F){
    data <- NULL
    if(type=="consensus"){
      predict.consensus <- predict(res, what="consensus")
      silhouette.consensus <- silhouette(res, what="consensus")
      # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
      # that matches the original sampleNames(res) order
      # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
      # order of the samples in consensusmap(res). It is just for displaying
      # Therefore, the merged data frame sampleNames(res) + a.predict.consensus is the final
      # consensus results.
      data <- data.frame(Sample_ID=sampleNames(res),
                         nmf_subtypes = predict.consensus,
                         sil_width = signif(silhouette.consensus[, "sil_width"], 3))
      # If we want to display as we see in consensusmap, we just need to reoder everything.
      # Now re-order data to match consensusmap sample order
      if(matchConseOrder){
        sample.order <- attributes(predict.consensus)$iOrd
        data <- data[sample.order, ]
      }
    }else if(type=="samples"){
      predict.samples <- predict(res, what="samples", prob=T)
      silhouette.samples <- silhouette(res, what="samples")
      data <- data.frame(Sample_ID=names(predict.samples$predict),
                         nmf_subtypes = predict.samples$predict,
                         sil_width = signif(silhouette.samples[, "sil_width"], 3),
                         prob = signif(predict.samples$prob, 3))
    }else{
      stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
    }
    return(data)
}
                            
write.table(nmf_extract_group(res2, type="consensus", matchConseOrder=T), "genes3clusters_consensusgroups.txt")
write.table(nmf_extract_group(res2, type="samples", matchConseOrder=T), "genes3clusters_samplegroups.txt")                           
                            
### obtain summary plots for these results
nmf_plot <- function(res, type="consensus", subsetRow=TRUE, save.image=F, hclustfun="average", silorder=F, add_original_name=T){
   if(save.image)
      pdf(res.pdf, width=18, height=15)

   if(type=="result"){
      print(plot(res))
   }else{
      si <- silhouette(res, what=type)

      if(type=="features"){
         if(silorder){
            basismap(res, Rowv = si, subsetRow=subsetRow)
         }else{
            basismap(res, subsetRow = subsetRow)
         }
      }else if(type=="samples"){
         if(silorder){
            coefmap(res, Colv = si)
         }else{
            coefmap(res)
         }
      }else if(type=="consensus"){
         if(add_original_name){
            colnames(res@consensus) <- sampleNames(res)
            rownames(res@consensus) <- sampleNames(res)
         }
         consensusmap(res, hclustfun=hclustfun)
      }
   }
   if(save.image)
      dev.off()
}

nmf_plot(res2, type="features", silorder=T, save.image=T)
nmf_plot(res2, type="consensus", silorder=T, save.image=T)
nmf_plot(res2, type="samples", silorder=T, save.image=T)
