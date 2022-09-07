#####======================================
#####1.Computation, Analysis and Prediction
#####======================================
#####Based on PRISM and CTRP2.0 drug sensitivity AUC value data, and CCLE expression profile data, potential therapeutic drugs in TCGA HNSCC subgroups were predicted. 
###install.packages("tidyverse")
library(tidyverse) # For reading files
###install.packages("ISOpureR")
library(ISOpureR) # For purification of expression profiles
###if (!requireNamespace("BiocManager", quietly = TRUE))
###install.packages("BiocManager")
###BiocManager::install("impute")
library(impute) # Use KNN to fill in drug susceptibility data
###install.packages("SimDesign")
library(SimDesign) # Information used to suppress the output of the susceptibility prediction process
library(ggplot2) # Drawing
library(cowplot) # Merge images
###Install the dependencies of the pRRophetic package
###BiocManager::install("sva")
###BiocManager::install("genefilter", force = TRUE)
###BiocManager::install("preprocessCore")
###install.packages("ridge")
###Use the Chinese domestic mirror installation package
###options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
###options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
###install.packages("pRRophetic_0.5.tar.gz", repos = NULL, type = "source")
###Use the downloaded local installation package
###install.packages("/home/tlq/anaconda3/Rpackages_by_tlq/pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
###install.packages("D:\\My_Win_home_tlq_data\\Rpackages_by_tlq\\pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
library(pRRophetic) # For drug susceptibility prediction

###===============================
###Test whether pRRophetic can be used normally
set.seed(1234)
###Load the data and draw a picture to see:Load the data and draw a picture to see:Load the data and draw a picture to see:Load the data and draw a picture to see:
data("bortezomibData")
pRRopheticQQplot("Bortezomib")
cvOut <- pRRopheticCV("Bortezomib", cvFold=5, testExprData=exprDataBortezomib)
###11683  gene identifiers overlap between the supplied expression matrices... 
###Found2batches
###Adjusting for0covariate(s) or covariate level(s)
###Standardizing Data across genes
###Fitting L/S model and finding priors
###Finding parametric adjustments
###Adjusting the Data
###1 of 5 iterations complete.
###2 of 5 iterations complete.
###3 of 5 iterations complete.
###4 of 5 iterations complete.
###5 of 5 iterations complete.
plot(cvOut)
summary(cvOut)
###Summary of cross-validation results:
###  Pearsons correlation: 0.4 , P =  4.45287272844977e-12 
###R-squared value: 0.16
###Estimated 95% confidence intervals: -4.23, 4.23
###Mean prediction error: 1.64
###With a model, you can make predictions:
predictedPtype <- pRRopheticPredict(exprDataBortezomib, "Bortezomib", selection=1)
###11683  gene identifiers overlap between the supplied expression matrices... 
###Found2batches
###Adjusting for0covariate(s) or covariate level(s)
###Standardizing Data across genes
###Fitting L/S model and finding priors
###Finding parametric adjustments
###Adjusting the Data
###2324 low variabilty genes filtered.
###Fitting Ridge Regression model... Done
###Calculating predicted phenotype...Done
###===============================


###Custom function to show progress
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}



#####=================================
####First, you need to put the following necessary common 6 files into the ../myresult_txt/ folder
##CTRP_AUC_raw.txt(CTRP_AUC.txt will be generated during operation),
##CTRP_ccl_anno.txt，
##CTRP_cpd_anno.txt,
##PRISM_AUC.txt,
##overlapTable27.txt
##CCLE_RNAseq_rsem_genes_tpm_20180929.txt
##Besides, you need to put the following necessary common 6 files into the ../myresult_txt/ folder
##HNSC_gene_expression_file02_tumor_no_col1name.tsv
##HNSC_gene_expression_file07t_mytargetgene_no_col1name.tsv

##HNSC_gene_expression_file02_tumor_no_col1name.tsv，Gene Expression Matrix. The myscore score and drug susceptibility prediction are calculated based on this expression matrix. Replace this expression matrix with your own data in practical applications.
##10gene.txt was used to calculate myscore score. 10-gene set in the model with the largest C-index value was considered as myscore.
##CTRP_AUC_raw.txt was used to make CTRP AUC matrix, raw data was from Correlating chemical sensitivity and basal gene expression reveals mechanism of action，2016 Nature Chemical Biology，Supplementary Data Set 1 2
##CTRP_ccl_anno.txt, CTRP_cpd_anno.txt, were respectively from Supplementary Data Set 1, Supplementary Data Set 2
##PRISM_AUC.txt was from https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
##CCLE_RNAseq_rsem_genes_tpm_20180929.txt was cell line expression profiles, as a training set for drug susceptibility prediction, which was from https://portals.broadinstitute.org/ccle/data.


### 1.Reading the expression profile of HNSCC
expr <- read.table("../myresult_txt/HNSC_gene_expression_file02_tumor_no_col1name.tsv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
normsam <- colnames(expr[,which(substr(colnames(expr),14,15) == "11")])
tumosam <- colnames(expr[,which(substr(colnames(expr),14,15) == "01")])
expr <- expr[-which(rownames(expr) %in% c("Bcr_patient_barcode", "Bcr_patient_typeid", "Bcr_patient_type")),]
### 2.Read profile(downloaded from cBioPortal)
Sinfo <- read.table("../myresult_txt/HNSC_gene_expression_file07t_mytargetgene_no_col1name.tsv", sep = "\t", header = T, row.names = 1)
Sinfo <- Sinfo[, which(names(Sinfo) %in% c("Cluster_subgroup", "Risk_score"))]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
###Standardize, process Risk_score between 0-1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Sinfo$Risk_score <- range01(Sinfo$Risk_score) 
Sinfo <- Sinfo[which(Sinfo$Cluster_subgroup %in% "High-risk"), ]
Sinfo$Cluster_subgroup <- gsub("High-risk", "1", Sinfo$Cluster_subgroup)
tumosam_High_risk <- rownames(Sinfo)
head(Sinfo)



###Expression profiles of tumor samples purified with ISOpureR
###Calculate the purified expression profile for correlation analysis with drug susceptibility data (this analysis requires the use of the original expression profile, logarithmic transformation is not possible)

###The running time is about a week, so it is not run here, and the original expression profile is used instead of downstream analysis

###If you want to run, change runpure <- F to runpure <- T

normexpr <- expr[,which(colnames(expr) %in% normsam)]
tumoexpr <- expr[,which(colnames(expr) %in% tumosam_High_risk)]
###normexpr <- expr[,normsam]
###tumoexpr <- expr[,tumosam]
###normexpr <- as.matrix(expr[,normsam])
###tumoexpr <- as.matrix(expr[,tumosam])

runpure <- F # Change this to T if you want to run
if(runpure) {
  set.seed(123)
  # Run ISOpureR Step 1 - Cancer Profile Estimation
  ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
  # For reproducible results, set the random seed
  set.seed(456);
  # Run ISOpureR Step 2 - Patient Profile Estimation
  ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
  pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
  pure.tumoexpr <- tumoexpr
}
write.table(pure.tumoexpr,"../myresult_txt/HNSC_gene_expression_file02_puretumoexpr_no_col1name.tsv",sep = "\t",row.names = T,col.names = T,quote = F)




###Drug susceptibility data preparation and preprocessing

### 1.Make CTRP AUC matrix and save it to CTRP_AUC.txt file
auc <- read.table("../myresult_txt/CTRP_AUC_raw.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 3
auc$comb <- paste(auc$master_cpd_id,auc$master_ccl_id,sep = "-")
auc <- apply(auc[,"area_under_curve",drop = F], 2, function(x) tapply(x, INDEX=factor(auc$comb), FUN=max, na.rm=TRUE)) # Duplicates take maximum AUC
auc <- as.data.frame(auc)
auc$master_cpd_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",1)
auc$master_ccl_id <- sapply(strsplit(rownames(auc),"-",fixed = T),"[",2)
auc <- reshape(auc, 
               direction = "wide",
               timevar = "master_cpd_id",
               idvar = "master_ccl_id")
colnames(auc) <- gsub("area_under_curve.","",colnames(auc),fixed = T)
ctrp.ccl.anno <- read.table("../myresult_txt/CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 1
ctrp.cpd.anno <- read.delim("../myresult_txt/CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # Supplementary Data Set 2
# Save to file
write.table(auc,"../myresult_txt/CTRP_AUC.txt",sep = "\t",row.names = F,col.names = T,quote = F)

### 2.Load the drug susceptibility AUC matrix and perform data processing
ctrp.auc <- read.table("../myresult_txt/CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim("../myresult_txt/PRISM_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # The data was from https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5] # The first 5 columns are cell line annotation information
prism.auc <- prism.auc[,-c(1:5)]

## a. Remove drugs with more than 20% missing values
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]

## b. Removed cell lines derived from haematopoietic_and_lymphoid_tissue in CTRP data
rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[which(ctrp.ccl.anno$ccle_primary_site == "haematopoietic_and_lymphoid_tissue"),"master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[setdiff(rownames(ctrp.auc),rmccl),]

## c. KNN fills missing values
ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data
## Warning in knnimp(x, k, maxmiss = rowmax, maxp = maxp): 30 rows with more than 50 % entries missing;
##  mean imputation used for these rows
prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data
## Warning in knnimp(x, k, maxmiss = rowmax, maxp = maxp): 1 rows with more than 50 % entries missing;
##  mean imputation used for these rows
## d. Data magnitude correction
ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) # Refer to Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))
### Drug susceptibility prediction
### Load expression profiles of CCLE cell lines as training set
ccl.expr <- read.table("../myresult_txt/CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

### Load gene annotation files for gene ID conversion
Ginfo <- read.table("../myresult_txt/overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

### Convert gene's ensembl ID to gene symbol
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]
###The following uses the calcPhenotype function in the pRRophetic package to estimate the drug response of each sample based on CTRP and PRISM, respectively.



###CTRP
###Data preparation

keepgene <- apply(ccl.expr, 1, mad) > 0.5 # Genes with valid expression values retained
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # reset cell line name
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() # Alternate cell line name
for (i in rownames(trainPtype)) {
  if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
    cat(i,"\n")
    ccl.miss <- c(ccl.miss, i) # No matching cell line
    ccl.name <- c(ccl.name, i) # Insert unmatched cell line
  } else {
    ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) # Insert matched cell line
  }
}
## CCL223 
## CCL1017 
## CCL1019 
## CCL1073 
## CCL1245 
## CCL241 
## CCL381 
## CCL382 
## CCL500 
## CCL895 
## CCL909 
## CCL990
cpd.name <- cpd.miss <- c() # alternate drug name
for (i in colnames(trainPtype)) {
  if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
    cat(i,"\n")
    cpd.miss <- c(cpd.miss, i) # No matching drug
    cpd.name <- c(cpd.name, i) # Insert unmatched drug
  } else {
    cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) # Insert matching drug
  }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # Remove unmatched cell lines
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # Remove unmatched drugs
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # Extraction of expressing and drug-sensitive cell lines
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

### Test dataset
pure.tumoexpr <- read.table("../myresult_txt/HNSC_gene_expression_file02_puretumoexpr_no_col1name.tsv",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5 # The purified test set was used to select genes with stable expression
###keepgene <- apply(as.numeric(unlist(pure.tumoexpr)), 1, mad) > 0.5  # The purified test set takes the gene with stable expression #If an error is reported, it can be solved by this method
testExpr <- log2(pure.tumoexpr[keepgene,] + 1) # logarithmic expression profile
### Take the genes common to the training set and the test set
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
###Ridge regression predicts drug sensitivity

outTab <- NULL
# The loop is slow, please be patient
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # Since the AUC of CTRP may have a value of 0, add a smaller value to prevent an error
  
  # Ridge regression predicts drug sensitivity
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),  ##"testExprData" must be a matrix.
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # antilog
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
## 5% 10% 15% 20% 25% 31% 36% 41% 46% 51% 56% 61% 66% 71% 76% 81% 86% 92% 97%
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
write.table(ctrp.pred.auc,file="../myresult_txt/HNSC_ctrp.pred.auc.tsv",sep="\t", quote=F, col.names=T, row.names=T)


###PRISM
###Data preparation

keepgene <- apply(ccl.expr, 1, mad) > 0.5
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
#colnames(trainPtype) <- sapply(strsplit(colnames(trainPtype)," (",fixed = T), "[",1)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

###Test dataset
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
###Ridge regression predicts drug sensitivity

outTab <- NULL
### The loop is slow, please be patient
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # Since the AUC of PRISM may have a value of 0, add a smaller value to prevent an error
  ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = as.matrix(testExpr),  ##"testExprData" must be a matrix.
                                  powerTransformPhenotype = F,
                                  selection = 1))
  ptypeOut <- 2^ptypeOut - 0.00001 # antilog
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
## 5% 10% 15% 20% 25% 31% 36% 41% 46% 51% 56% 61% 66% 71% 76% 81% 86% 92% 97%
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
prism.pred.auc <- outTab
write.table(prism.pred.auc,file="../myresult_txt/HNSC_prism.pred.auc.tsv",sep="\t", quote=F, col.names=T, row.names=T)

###Identify potential drug targets
top.myscore <- Sinfo[Sinfo$Risk_score >= quantile(Sinfo$Risk_score,probs = seq(0,1,0.1))[10],] # Samples that define the upper decile
bot.myscore <- Sinfo[Sinfo$Risk_score <= quantile(Sinfo$Risk_score,probs = seq(0,1,0.1))[2],] # Define a sample of the lower decile

###1.Differential susceptibility analysis
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.myscore)])) # AUC mean of upper decile
  b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.myscore)])) # AUC mean of lower decile
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
## 5% 10% 15% 20% 25% 31% 36% 41% 46% 51% 56% 61% 66% 71% 76% 81% 86% 92% 97%
candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.1] # Here I adjusted the threshold to control the number of results
write.table(candidate.ctrp,file="../myresult_txt/HNSC_candidate.ctrp.tsv",sep="\t", quote=F, col.names=T, row.names=T)

prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,rownames(top.myscore)])) # AUC mean of upper decile
  b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.myscore)])) # AUC mean of lower decile
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
## 5% 10% 15% 20% 25% 30% 35% 40% 46% 51% 56% 61% 66% 71% 76% 81% 86% 91% 96%
candidate.prism <- prism.log2fc[prism.log2fc > 0.1] # Here I adjusted the threshold to control the number of results
write.table(candidate.prism,file="../myresult_txt/HNSC_candidate.prism.tsv",sep="\t", quote=F, col.names=T, row.names=T)

###2.Spearman correlation analysis, used to draw the left graph
ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$Risk_score)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
## 5% 10% 15% 20% 25% 31% 36% 41% 46% 51% 56% 61% 66% 71% 76% 81% 86% 92% 97%
write.table(ctrp.cor,file="../myresult_txt/HNSC_candidate.ctrp.cor.tsv",sep="\t", quote=F, col.names=T, row.names=T)
###Only take negative correlations
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.3]  # Here I adjusted the threshold to control the number of results
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))
###Take both positive and negative correlations
candidate.ctrp3 <- ctrp.cor[(ctrp.cor < 0.3) | (ctrp.cor > 0.3)]  # Here I adjusted the threshold to control the number of results
ctrp.candidate2 <- intersect(names(candidate.ctrp),names(candidate.ctrp3))


prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$Risk_score)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
## 5% 10% 15% 20% 25% 30% 35% 40% 46% 51% 56% 61% 66% 71% 76% 81% 86% 91% 96%
write.table(prism.cor,file="../myresult_txt/HNSC_candidate.prism.cor.tsv",sep="\t", quote=F, col.names=T, row.names=T)
###Only take negative correlations
candidate.prism2 <- prism.cor[prism.cor < -0.3]  # Here I adjusted the threshold to control the number of results
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
###Take both positive and negative correlations
candidate.prism3 <- prism.cor[(prism.cor < -0.3) | (prism.cor > 0.3)]  # Here I adjusted the threshold to control the number of results
prism.candidate2 <- intersect(names(candidate.prism),names(candidate.prism3))

###Start drawing (only take negative correlation)

###1. Correlation Plot
# Set color
darkblue <- "#0772B9"
lightblue <- "#48C8EF"

cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #left and right blank
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #left and right blank
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

###2.Box plot
ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.myscore)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.myscore)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High Score","Low Score"),c(nrow(top.myscore),nrow(bot.myscore))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
                label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.myscore)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(bot.myscore)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High Score","Low Score"),c(nrow(top.myscore),nrow(bot.myscore))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
      legend.position = "bottom",
      legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

###3. Merge images
plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2)) #The width ratio of the left and right columns


ggsave(filename = "../myresult_pdf/HNSC_drug target_negative.pdf",width = 8,height = 8)




###Start drawing (take both positive and negative correlations)
pdf("../myresult_pdf/HNSC_drug target.pdf",width = 3.3,height = 3.3)
###1. Correlation Plot
# set color
darkblue <- "#E64B35FF"
lightblue <- "#4DBBD5FF"
###darkblue <- "#0772B9"
###lightblue <- "#48C8EF"

cor.data <- data.frame(drug = ctrp.candidate2,
                       r = ctrp.cor[ctrp.candidate2],
                       p = -log10(ctrp.cor.p[ctrp.candidate2]))

cor.data <- data.frame(drug = c("paclitaxel", "doxorubicin", "gemcitabine", "vincristine", "SB-743921", "clofarabine", "rigosertib"),
                      r = c(-0.4102130, -0.3879626, -0.3688861, -0.3532887, -0.3454607, -0.3316535, -0.3095948),
                       p = c(6, 6, 5, 5, 5, 4, 4))
###                  drug          r        p
###paclitaxel   paclitaxel -0.4202130 5.553666
###doxorubicin doxorubicin -0.3879626 4.961512
###gemcitabine gemcitabine -0.3688861 3.215833
###vincristine vincristine -0.3532887 3.053565
###SB-743921     SB-743921 -0.3454607 2.702371
###clofarabine clofarabine -0.3316535 2.855352
###rigosertib   rigosertib -0.3095948 2.393636
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,6)) +
  scale_x_reverse(breaks = c(0, -0.1, -0.2, -0.3, -0.4, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #left and right blank
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "top", 
        axis.line.y = element_blank())
print(p1)

cor.data <- data.frame(drug = prism.candidate2,
                       r = prism.cor[prism.candidate2],
                       p = -log10(prism.cor.p[prism.candidate2]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,6)) +
  scale_x_reverse(breaks = c(0, -0.1, -0.2, -0.3, -0.4, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #left and right blank
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "top", 
        axis.line.y = element_blank())
print(p2)

###2.Box plot
ctrp.boxdata <- NULL
for (d in ctrp.candidate2) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.myscore)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.myscore)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High risk score","Low risk score"),c(nrow(top.myscore),nrow(bot.myscore))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
                label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 10),
        legend.position = "top",
        legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="black", inherit.aes = F)
print(p3)

prism.boxdata <- NULL
for (d in prism.candidate2) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.myscore)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(bot.myscore)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High risk score","Low risk score"),c(nrow(top.myscore),nrow(bot.myscore))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 10),
      legend.position = "top",
      legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="black", inherit.aes = F)
print(p4)

dev.off()