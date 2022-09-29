#2nd Analysis/2nd code

BiocManager::install('org.Hs.eg.db') # Annotation of the human genome
library(org.Hs.eg.db)
BiocManager::install("AnnotationDbi")
BiocManager::install("hugene20sttranscriptcluster.db")
library(hugene20sttranscriptcluster.db)

install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("pheatmap")


library(tidyverse)
library(dplyr)
library(GEOquery)
library(limma) #linear models and differential expression for microarray data.
library(oligo)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


gse56481 <- getGEO('GSE56481',getGPL= TRUE, GSEMatrix=TRUE)[[1]]
gse56481

names(pData(gse56481))[38] <- 'diagnosis'
names(pData(gse56481))[39] <- 'facs_cell_type'
#names(pData(gse56481))

pd <- pData(gse56481)


pd$diagnosis <- as.factor(pd$diagnosis)
levels(pd$diagnosis) <- c("GPA","Control")

#pd$diagnosis

pd$group <- as.factor(paste(pd$diagnosis,pd$facs_cell_type))
pd$group


levels(pd$group) <- c ("Control.CD4positive", "Control.CD4_CD8positive", "Control.CD8positive",  "GPA.CD4positive" ,  "GPA.CD4_CD8positive" , "GPA.CD8positive" )

design <- model.matrix(~ 0 + pd$group)
colnames(design) <- levels(pd$group)
design

contrasts_matrix <- makeContrasts(de_CD4_CD8positive = GPA.CD4_CD8positive - Control.CD4_CD8positive ,
                                  de_CD4positive = GPA.CD4positive - Control.CD4positive ,
                                  de_CD8positive = GPA.CD8positive - Control.CD8positive , 
                                  levels = design)
 
#dim(gse56481)   
#dim(design) 
                                  
gse56481_fit <- lmFit(gse56481,design)# fits a linear model using a gene expression object, fits a linear model for every single gene
gse56481_fit2 <- contrasts.fit(gse56481_fit,contrasts=contrasts_matrix)
gse56481_fit2 <- eBayes(gse56481_fit2)# modifying our estimates by a prior expectation

summary(decideTests(gse56481_fit2, p.value=0.05,lfc=1))
head(decideTests(gse56481_fit2,lfc=1),10)

topTable(gse56481_fit2)


#ps <- rownames(topTable(gse56481_fit2))

#ls('package:hugene20sttranscriptcluster.db')
#columns(hugene20sttranscriptcluster.db)
#keytypes(hugene20sttranscriptcluster.db)
#head(keys(hugene20sttranscriptcluster.db,keytype="PROBEID"))
#AnnotationDbi::select(hugene20sttranscriptcluster.db, ps ,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

ps2 <- topTable(gse56481_fit2,number=Inf,p.value = 0.05,lfc=2)

ps2_up <- rownames(ps2[ps2$de_CD4_CD8positive > 0,])
ps2_up
df_up <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df_up,GENENAME=stringr::str_trunc(GENENAME,30))

ps2_down <- rownames(ps2[ps2$de_CD4_CD8positive < 0,])
ps2_down
df_down <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_down,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
dplyr::mutate(df_down,GENENAME=stringr::str_trunc(GENENAME,30))

         

interesting_genes <- topTable(gse56481_fit2,number=Inf,p.value = 0.05,lfc=2)
interesting_genes

#volcano plot
volcanoplot(gse56481_fit2, coef=2, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))

points(interesting_genes[['de_CD4_CD8positive']],-log10(interesting_genes[['P.Value']]),col='red') #differentially expressed transcripts of CD4+CD8+ cells, GPA.CD4_CD8positive - Control.CD4_CD8positive
points(interesting_genes[['de_CD4positive']],-log10(interesting_genes[['P.Value']]),col='blue') #DEGs of CD4+ cells
points(interesting_genes[['de_CD8positive']],-log10(interesting_genes[['P.Value']]),col='green') #DEGs of CD8+ cells


#Heatmap
eset_of_interest <- gse56481[rownames(interesting_genes),]

pheatmap(exprs(eset_of_interest), cellwidth = 13, cellheight = 3.2, fontsize_row = 4 , labels_col=pd$group )

















































