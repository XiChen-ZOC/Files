rm(list = ls())
exprSet <- data.table::fread("data/KC2_rawcount.csv",data.table = F)

library(tidyverse)
exprSet <- exprSet %>%
  distinct(gene,.keep_all = T)
  
metadata <- data.frame(row.names = names(exprSet)[-1],
                       'sample' = names(exprSet)[-1],
                       'group' = c(rep('NC',3), rep('KC',12)))

save(metadata,file = "output/metadata_KC2.Rdata")

library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet,
                             colData=metadata,
                             design=~group,
                             tidy=TRUE)
nrow(dds)
rownames(dds)

# counts(dds)
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "group", returnData=F)

library(clusterProfiler)
library(org.Hs.eg.db)
gene <- 'xxxx'
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL"),
                OrgDb = org.Hs.eg.db);gene.df

plotCounts(dds, gene = gene.df$SYMBOL, intgroup=c("group"), returnData=T)
plotCounts(dds, gene = gene.df$SYMBOL, intgroup=c("group"), returnData=F)

exprSet_vst <- as.data.frame(assay(vsd))
# test <- exprSet_vst[1:10,1:10]
save(exprSet_vst,file = "output/exprSet_vst_KC2.Rdata")

dds <- DESeq(dds)

contrast=c("group", "KC", "NC")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
# install.packages('ashr')
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1, type="ashr")
plotMA(dd2, ylim=c(-5,5))

library(dplyr)
library(tibble)
library(tidyr)
res <- dd2 %>%
  as.data.frame() %>%
  rownames_to_column("gene")

library(AnnotationDbi)
library(org.Hs.eg.db)
res$gene_id <- mapIds(org.Hs.eg.db,
                      keys=res$gene,
                      column="ENSEMBL",
                      keytype="SYMBOL",
                      multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

colnames(res) <- c("gene","baseMean","logFC","lfcSE","P.Value","adj.P.Val","gene_id","entrez")

save(res,file = "output/DEseq2_Diff_KC2.Rdata")
