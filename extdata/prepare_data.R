library(slingshot)
library(clusterExperiment)

load("extdata/regenK5_1_se_filtqc_idfiltyes_20190530_160100.Rda")

oe <- readRDS("extdata/slingshot_OE.rds")
load("extdata/regenK5_scone_none,fq,ruv_k=1,no_bio,batch_rsec_adjP_mergecutoff_0.01_20190609_085359.Rda")
idx <- readRDS("extdata/subset_index.rds")

cl3 <- cl2[,idx]
rm(cl2)
colData(cl3)

cellidx <- unlist(apply(oe@clusterLabels[,oe@lineages$Lineage1], 2, function(x) which(x==1)))

cl4 <- cl3[,cellidx]

# tfs <- read.table("https://raw.githubusercontent.com/rufletch/p63-HBC-diff/master/ref/ATFDB_mm_TF.txt",
#                   stringsAsFactors = FALSE)[,1]
#
# genes <- intersect(rownames(cl4), tfs)
# length(genes)
# head(genes)
#
# cl5 <- cl4[genes,]
se <- se_filtered[,colData(cl4)$samples]

datamat <- t(assay(se))
any(rowSums(datamat)==0)
any(colSums(datamat)==0)

write.table(datamat, file = "extdata/oe_data_10x_neuronal.txt", sep="\t", quote = FALSE, row.names = TRUE)

tt <- droplevels(se$expt)
tt <- model.matrix(~tt - 1)
write.table(tt, file = "extdata/oe_data_10x_neuronal_covariates.txt", sep="\t", quote = FALSE, row.names = FALSE)
