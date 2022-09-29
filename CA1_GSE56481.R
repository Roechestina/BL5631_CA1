library(dplyr)
library(GEOquery)
library(hugene20sttranscriptcluster.db)
library(limma)
library(oligo)
library(oligoClasses)
library(RColorBrewer)
library(tidyverse)


#methods(class=class(gse56481))

gse56481 <- getGEO('GSE56481')
gse56481 <- gse56481[[1]]

# Guarantee that CEL files are in the correct order as the experimental data 
# we have from getGEO
pd <- pData(gse56481)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail, 1)
# Read the cel files
gse56481_celdata <- read.celfiles(paste0('GSE56481_RAW/',pd$cel_file),
                                  phenoData=phenoData(gse56481))

# Look for variables of interest
varLabels(gse56481) 
# Extract data of interest
pData(gse56481_celdata)[,c('geo_accession', 'cell type:ch1', 'diagnosis:ch1')]

# Data processing with RMA
gse56481_eset <- rma(gse56481_celdata)

# Identify differentially expressed genes
design <- model.matrix(~ gse56481_eset[['diagnosis:ch1']])
colnames(design)[2] <- 'Control'

fit <- lmFit(gse56481_eset, design)
fitted.ebayes <- eBayes(fit)
topTable(fitted.ebayes)
summary(decideTests(fitted.ebayes[, "Control"], lfc=1))

# Contrast matrix
design <- model.matrix( ~ 0 + gse56481_eset[['diagnosis:ch1']])
colnames(design) <- c('GPA', "Control")
contrast_matrix <- makeContrasts(GPA - Control, Control - GPA, levels=design)
contrast_matrix

fit <- lmFit(gse56481_eset, design)
fit2 <- contrasts.fit(fit, contrasts=contrast_matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2, lfc=1))

# Annotation of genomic data 
# using hugene20sttranscriptcluster.db (pd.hugene.2.0.st)
ps <- rownames(topTable(fitted.ebayes))

ls('package:hugene20sttranscriptcluster.db')
unlist(mget(ps,hugene20sttranscriptclusterSYMBOL))

columns(hugene20sttranscriptcluster.db)
keytypes(hugene20sttranscriptcluster.db)

head(keys(hugene20sttranscriptcluster.db, keytype="PROBEID"))
AnnotationDbi::select(hugene20sttranscriptcluster.db, ps, 
                      c("SYMBOL", "ENTREZID", "GENENAME"), keytype="PROBEID")

# retrieve all genes that are differentially expressed 
# with a adjusted p-value of less than 0.05 and log fold change one
ps2 <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1)
ps2_up <- rownames(ps2[ps2$logFC > 0, ]) # up-regulated genes
ps2_down <- rownames(ps2[ps2$logFC < 0, ]) # down-regulated genes
df_up <- AnnotationDbi::select(hugene20sttranscriptcluster.db, ps2_up,  
                            c("SYMBOL", "ENTREZID", "GENENAME"), keytype="PROBEID")
df_down <- AnnotationDbi::select(hugene20sttranscriptcluster.db, ps2_down, 
                                 c("SYMBOL", "ENTREZID", "GENENAME"), keytype="PROBEID")
dplyr::mutate(df_up,GENENAME=stringr::str_trunc(GENENAME,30))
dplyr::mutate(df_down,GENENAME=stringr::str_trunc(GENENAME,30))

# Data visualization
interesting_genes <- topTable(fitted.ebayes, number=Inf, p.value = 0.05, lfc=2)
volcanoplot(fitted.ebayes, coef=2, main=sprintf("%d features pass our cutoffs", 
                                                nrow(interesting_genes)))
points(interesting_genes[['logFC']], -log10(interesting_genes[['P.Value']]), col='red')

eset_of_interest <- gse56481_eset[rownames(interesting_genes), ]

heatmap(exprs(eset_of_interest), 
        labCol=gse56481_eset[['diagnosis:ch1']] , labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")), 
        distfun   = function(x) as.dist(1-cor(t(x))), revC = TRUE)
