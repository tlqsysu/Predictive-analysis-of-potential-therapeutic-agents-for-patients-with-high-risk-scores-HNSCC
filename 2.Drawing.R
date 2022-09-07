#####======================================
#####2.Drawing
#####======================================
#####We used data from multiple sources to corroborate each other to illustrate the potential value of the screened compounds/drugs. We designed this diagram to neatly display evidences 1, 2, 3, and 4 at the same time, while also seeing the overlap.

###Evidence I and II：We found 16 candidate compounds identified showed a higher drug sensitivity in high-risk HNSCC patients through PRISM and CTRP.
###Evidence III：We first used the CMap analysis to find compounds of which gene expression patterns were oppositional to the HNSCC-specific expression patterns (i.e. gene expression increased in tumor tissues but decreased by treatment of certain compounds).
###Evidence IV：Secondly, fold-change differences of the expression levels (including mRNA- and protein-level) of candidates’ drug targets between tumor and normal tissue were calculated, and a higher fold change value indicated a greater potential of candidate agent for HNSCC treatment.
###Evidence V：Thirdly, a comprehensive literature search was performed in PubMed (https://www.ncbi.nlm.nih. gov/pubmed/) to find out the experimental and clinical evidence of candidate compounds in treating HNSCC.

###The following will take you to implement CMap analysis and drawing.

library(devtools)
### Loading required package: usethis
library(ComplexHeatmap) # for drawing heatmaps
### Loading required package: grid
### ========================================
### ComplexHeatmap version 2.4.3
### Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
### Github page: https://github.com/jokergoo/ComplexHeatmap
### Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
### 
### If you use it in published research, please cite:
### Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
###   genomic data. Bioinformatics 2016.
### 
### This message can be suppressed by:
###   suppressPackageStartupMessages(library(ComplexHeatmap))
### ========================================
library(circlize) # for color matching
### ========================================
### circlize version 0.4.11
### CRAN page: https://cran.r-project.org/package=circlize
### Github page: https://github.com/jokergoo/circlize
### Documentation: https://jokergoo.github.io/circlize_book/book/
### 
### If you use it in published research, please cite:
### Gu, Z. circlize implements and enhances circular visualization
###   in R. Bioinformatics 2014.
### 
### This message can be suppressed by:
###   suppressPackageStartupMessages(library(circlize))
### ========================================
library(tidyverse) # for reading various files
### ─ Attaching packages ──────────────────── tidyverse 1.3.0 ─
### ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
### ✓ tibble  3.0.4     ✓ dplyr   1.0.2
### ✓ tidyr   1.1.2     ✓ stringr 1.4.0
### ✓ readr   1.4.0     ✓ forcats 0.5.0
### ─ Conflicts ───────────────────── tidyverse_conflicts() ─
### x dplyr::filter() masks stats::filter()
### x dplyr::lag()    masks stats::lag()
library(limma) # for differential expression

Sys.setenv(LANGUAGE = "en") #Display English error message
options(stringsAsFactors = FALSE) #Forbid chr to be converted into factor

####=================================
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

delete_col1name("../myresult_txt/HNSC_gene_expression_file02.tsv", "Hybridization_REF", "../myresult_txt/HNSC_gene_expression_file02_no_col1name.tsv")

mydata <- read.table("../myresult_txt/HNSC_gene_expression_file02_no_col1name.tsv",sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
mydata <- mydata[-which(rownames(mydata) %in% c("Bcr_patient_barcode", "Bcr_patient_typeid", "Bcr_patient_type")),]
write.table(mydata,file="../myresult_txt/HNSC_gene_expression_file02_no_col1name_expr.tsv",sep="\t", quote=F, col.names=T, row.names=T)

### 1.Reading the expression profile of HNSCC
expr <- read.table("../myresult_txt/HNSC_gene_expression_file02_no_col1name_expr.tsv",sep = "\t", header = T, row.names = 1,check.names = F,stringsAsFactors = F)
normsam <- colnames(expr[,which(substr(colnames(expr),14,15) == "11")])
tumosam <- colnames(expr[,which(substr(colnames(expr),14,15) == "01")])

### 2.Read profile(downloaded from cBioPortal)
Sinfo <- read.table("../myresult_txt/HNSC_gene_expression_file07t_mytargetgene.tsv", sep = "\t", header = T, row.names = 1, check.names = F,stringsAsFactors = F)
Sinfo <- Sinfo[, which(names(Sinfo) %in% c("Cluster_subgroup", "Risk_score"))]
###Standardize, process Risk_score between 0-1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Sinfo$Risk_score <- range01(Sinfo$Risk_score) 
Sinfo <- Sinfo[which(Sinfo$Cluster_subgroup %in% "High-risk"), ]
Sinfo$Cluster_subgroup <- gsub("High-risk", "1", Sinfo$Cluster_subgroup)
tumosam_High_risk <- rownames(Sinfo)
head(Sinfo)

pd <- data.frame(Samples = c(tumosam_High_risk,normsam),
                 Group = rep(c("tumor","normal"),c(length(tumosam_High_risk),length(normsam))),
                 stringsAsFactors = FALSE)
design <-model.matrix(~ -1 + factor(pd$Group, levels = c("tumor","normal")))
colnames(design) <- c("tumor","normal")

gset <- log2(expr[,pd$Samples] + 1)
fit <- limma::lmFit(gset, design = design);
contrastsMatrix <- limma::makeContrasts(tumor - normal, levels = c("tumor", "normal"))
fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
fit2 <- limma::eBayes(fit2, 0.01)
## Warning: Zero sample variances detected, have been offset away from zero
resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
resData$id <- rownames(resData)
colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
resData$fc <- 2^resData$log2fc
resData <- resData[order(resData$padj),c("id","fc","log2fc","pvalue","padj")]

# Save to file
write.table(resData, file = "../myresult_txt/HNSCC_High_risk_output_degs.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

###CMap Analysis
###A lot of genes with the most significant fold changes (some up-regulated genes and some down-regulated genes) were submitted to the CMap website https://clue.io/query

ngene <- 150
degs <- na.omit(resData[order(resData$log2fc,decreasing = T),])
updegs <- rownames(degs)[1:ngene]
dndegs <- rownames(degs)[(nrow(degs)-ngene + 1):nrow(degs)]
cmap.input <- data.frame(up = updegs,
                         dn = dndegs,
                         stringsAsFactors = F)
write.table(cmap.input,"../myresult_txt/HNSCC_High_risk_CMap_input.tsv",sep = "\t",row.names = F,col.names = T,quote = F)
cmap.input.up <- data.frame(up = updegs,
                         stringsAsFactors = F)
write.table(cmap.input.up,"../myresult_txt/HNSCC_High_risk_CMap_input_upgenes.txt",sep = "\t",row.names = F,col.names = T,quote = F)
cmap.input.dn <- data.frame(dn = dndegs,
                         stringsAsFactors = F)
write.table(cmap.input.dn,"../myresult_txt/HNSCC_High_risk_CMap_input_dngenes.txt",sep = "\t",row.names = F,col.names = T,quote = F)



###CMap website usage
###Based on HNSCC_High_risk_CMap_input.tsv, by the following method, HNSCC_High_risk_CMap_export.txt was obtained and used for the next drawing



###start drawing
###The figure is divided into three parts: left, middle and right. The left block is equivalent to the annotation, the middle is the gene change fold, and the right is the CMap score.

###Fill in the annotation, log2FC, and CMap score calculated earlier.

###Here, the value is directly written into the matrix matrix, and in actual use, it can be assigned by reading the file storing the corresponding value.

cmap.res <- read.delim("../myresult_txt/HNSCC_High_risk_CMap_export.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

###The situation of CTRP database calculation and drawing
mydrug_CTRP <- c("paclitaxel", "doxorubicin", "gemcitabine", "vincristine", "SB-743921", "clofarabine", "rigosertib") 
mydrug <- mydrug_CTRP
selres <- cmap.res[which(cmap.res$Name %in% intersect(cmap.res$Name,mydrug)),]
print(selres) 

###The following is to build a matrix based on this result for drawing a heat map, and prepare separately:

###dummy data [for mapping colors],dt
###display label,lb
###The label corresponds to the color,cl
###left block
###background color
dt1 <- matrix(c(1,1, # 1 is dark, 0 is light
                1,1,
                1,1,
                1,1,
                0,0,
                0,0,
                1,1),
              ncol = 2,
              byrow = T,
              dimnames = list(mydrug, c("Clinical status","Experimental evidence")))

# text label
lb1 <- matrix(c("Phase3\nCompleted","Present",
                "Phase2\nCompleted","Present",
                "Phase2\nCompleted","Present",
                "Phase3\nCompleted","Present",
                "Preclinical","Absent",
                "Preclinical","Absent",
                "Phase2\nCompleted","Present"),
              ncol = 2,
              byrow = T)

# text color
cl1 <- matrix(c("white","white", #1 is dark background, corresponding to white font, 0 is light background, corresponding to black font
                "white","white",
                "white","white",
                "white","white",
                "black","black",
                "black","black",
                "white","white"),
              ncol = 2,
              byrow = T)

# drawing
hm1 <- pheatmap(mat = dt1,
                color = c("#B7E4ED", "#019AC9"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb1,
                number_color = cl1,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

# middle block
# background color
dt2 <- matrix(c(-1,-1, # 1 is the darkest color, 0 is the middle color (white) (NA), -1 is the lightest color
                1,1,
                -1,-1,
                1,-1,
                0,0,
                -1,-1,
                1,1),
              ncol = 2,
              byrow = T,
              dimnames = list(mydrug, c("Log2FC of mRNA expression","Log2FC of protein expression")))

# text label
lb2 <- matrix(c("0.09","0.04",
                "1.88","0.65",
                "0.19","0.07",
                "1.22","0.43",
                "NA","NA",
                "0.56","0.34",
                "2.03","0.87"),
              ncol = 2,
              byrow = T)
# text color
cl2 <- matrix(c("black","black", #1 is dark background, corresponding to white font, 0,-1 is light background, corresponding to black font
                "white","white",
                "black","black",
                "white","black",
                "black","black",
                "black","black",
                "white","white"),
              ncol = 2,
              byrow = T)
# drawing
hm2 <- pheatmap(mat = dt2,
                color = c( myblue, "#EFFAE8", myred ), #corresponds to -1,0,1
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb2,
                number_color = cl2,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

# right block
# background color
dt3 <- matrix(c(-1,  #1 is the darkest color, 0 is the middle color (white) (NA), -1 is the lightest color
                1,  
                -1,
                -1, 
                0,
                1,
                0),
              ncol = 1,
              byrow = T,
              dimnames = list(mydrug, c("CMap score")))

# text label
# Fill in according to the result of selres
# NA because there are two drugs that do not have CMap results
lb3 <- matrix(c("35.39",
                "-97.39",
                "-62",
                "54.47",
                "NA",
                "-98.03",
                "NA"),
              ncol = 1,
              byrow = T)

# text color
cl3 <- matrix(c("black",
                "white",
                "black",
                "black",
                "black",
                "white",
                "black"),
              ncol = 1,
              byrow = T)

# drawing
hm3 <- pheatmap(mat = dt3,
                color = c("#FDCEB9", "#EFFAE8", "#F26E5F"),  ###对应着-1,0,1###深红色###白色###浅红色
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb3,
                number_color = cl3,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

## Horizontal merge heatmap
hm <- hm1 + hm2 + hm3

pdf(file = "../myresult_pdf/HNSCC_High_risk_heatmap_CTRP.pdf",width = 8.0,height = 8.0)
draw(hm)
dev.off()
###Note that the entire picture is drawn by heat map stitching AI. The pdf file output here is a vector diagram. Please use AI to draw the remaining details and legends of the picture.



###The situation of PRISM database calculation and drawing
mydrug_PRISM <- c("paclitaxel", "docetaxel", "NVP-AUY922", "dasatinib", "epothilone-b", "talazoparib", "topotecan", "rubitecan", "vinblastine") 
mydrug <- mydrug_PRISM
selres <- cmap.res[which(cmap.res$Name %in% intersect(cmap.res$Name,mydrug)),]
print(selres) 

###he following is to build a matrix based on this result for drawing a heat map, and prepare separately:

###dummy data [for mapping colors],dt
###display label,lb
###The label corresponds to the color,cl
###left block
###background color
dt1 <- matrix(c(1,1, # 1 is dark, 0 is light
                1,1,
                0,0,
                1,1,
                1,1,
                1,1,
                1,1,
                0,0,
                0,0),
              ncol = 2,
              byrow = T,
              dimnames = list(mydrug, c("Clinical status","Experimental evidence")))

# text label
lb1 <- matrix(c("Phase3\nCompleted","Present",
                "Phase3\nCompleted","Present",
                "Preclinical","Absent",
                "Phase1\nCompleted","Present",
                "Phase1\nCompleted","Present",
                "Phase2\nCompleted","Present",
                "Phase2\nCompleted","Present",
                "Preclinical","Absent",
                "Preclinical","Absent"),
              ncol = 2,
              byrow = T)

# text color
cl1 <- matrix(c("white","white",
                "white","white",
                "black","black",
                "white","white",
                "white","white",
                "white","white",
                "white","white",
                "black","black",
                "black","black"),
              ncol = 2,
              byrow = T)

# drawing
hm1 <- pheatmap(mat = dt1,
                color = c("#B7E4ED", "#019AC9"),
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb1,
                number_color = cl1,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)


# middle block
# background color
dt2 <- matrix(c(-1,-1,  # 1 is the darkest color, 0 is the middle color (white) (NA), -1 is the lightest color
                -1,-1,
                0,0,
                -1,-1,
                1,1,
                -1,-1,
                1,1,
                -1,-1,
                1,1),
              ncol = 2,
              byrow = T,
              dimnames = list(mydrug, c("Log2FC of mRNA expression","Log2FC of protein expression")))

# text label
lb2 <- matrix(c("0.09","0.03",
                "0.09","0.04",
                "NA","NA",
                "0.97","0.47",
                "1.22","0.51",
                "0.86","0.44",
                "1.1","0.53",
                "0.36","0.17",
                "1.22","0.67"),
              ncol = 2,
              byrow = T)
#  text color
cl2 <- matrix(c("black","black",
                "black","black",
                "black","black",
                "black","black",
                "white","white",
                "black","black",
                "white","white",
                "black","black",
                "white","white"),
              ncol = 2,
              byrow = T)
# drawing
hm2 <- pheatmap(mat = dt2,
                color = c(myblue, "#EFFAE8", myred), ###corresponds to -1,0,1
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb2,
                number_color = cl2,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

# right block
# background color
dt3 <- matrix(c(-1, # 1 is the darkest color, 0 is the middle color (white) (NA), -1 is the lightest color
                -1, 
                1,
                -1, 
                0,
                0,
                1,
                0,
                -1),
              ncol = 1,
              byrow = T,
              dimnames = list(mydrug, c("CMap score")))

# text label
# Fill in according to the result of selres
# NA because there are two drugs that do not have CMap results
lb3 <- matrix(c("35.39",
                "76.26",
                "-96.12",
                "-60.24",
                "NA",
                "NA",
                "-98.17",
                 "NA",
                "11.27"),
              ncol = 1,
              byrow = T)

# text color
cl3 <- matrix(c("black",
                "black",
                "white",
                "black",
                "black",
                "black",
                "white",
                "black",
                "black"),
              ncol = 1,
              byrow = T)

# drawing
hm3 <- pheatmap(mat = dt3,
                color = c("#FDCEB9", "#EFFAE8", "#F26E5F"), ###corresponds to -1,0,1###dark red ###white ###light red
                cluster_cols = F,
                cluster_rows = F,
                show_rownames = T,
                show_colnames = T,
                display_numbers = lb3,
                number_color = cl3,
                fontsize_number = 10,
                border_color = "white",
                cellwidth = 50,
                cellheight = 50,
                legend = F)

## Horizontal merge heatmap
hm <- hm1 + hm2 + hm3 

pdf(file = "../myresult_pdf/HNSCC_High_risk_heatmap_PRISM.pdf",width = 8.0,height = 8.0)
draw(hm)
dev.off()
###Note that the entire picture is drawn by heat map stitching AI. The pdf file output here is a vector diagram. Please use AI to draw the remaining details and legends of the picture.