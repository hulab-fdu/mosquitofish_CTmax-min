library(DESeq2)

############DEG analysis using MAX as an example########

#Read in sample table and Trim counts table
counts_star=read.csv("htseq_count_table_1.csv", stringsAsFactors = F)
rownames(counts_star)=counts_star$gene
counts_star=counts_star[,-1]
noint_star = rownames(counts_star) %in% c("__no_feature","__ambiguous","__too_low_aQual", "__not_aligned","__alignment_not_unique")
samples_star=read.csv("samples_ren.csv", stringsAsFactors = F)
samples_star=samples_star[samples_star$libraryname %in% colnames(counts_star),]

# remove non-informative rows and weakly expressed genes in all samples
counts_star=counts_star[!noint_star,]
colData=data.frame(samples_star[,c("group")])
rownames(colData)=samples_star$libraryname 
colnames(colData)="group"
colnames(counts_star)==rownames(colData)

dds=DESeqDataSetFromMatrix(countData = counts_star, 
                           colData = colData,
                           design = ~ group) # design here is just fake to satisfy function

keep_gene=rowSums(counts(dds)>=1)>=nrow(colData) # remove genes with less than 1 read in all samples
dds=dds[keep_gene,]

#Trim sample table for all experimental pairs
keep_samples_max=grep("MAX", samples_star$condition)
samples_max=samples_star[keep_samples_max,]

keep_samples_min=grep("MIN", samples_star$condition)
samples_min=samples_star[keep_samples_min,]

#Trim counts table for DE analysis
max_keep=grep("max", colnames(counts_star)) 
max=dds[,max_keep]

min_keep=grep("min",colnames(counts_star))
min=dds[,min_keep]

# DE analysis of CTmax vs Control
colData_max=data.frame(samples_max[,c("group")])
rownames(colData_max)=samples_max$libraryname
colnames(colData_max)="group"

colnames(max)==rownames(colData_max) # check if the count table has the same order as sample table

max$group = factor(max$group, levels = c("Max_C","Max_T"))

# Test for DE genes
ddsMAX=DESeq(max)
resMAX=results(ddsMAX, contrast = c("group", "Max_T", "Max_C"))
head(resMAX)
exp_level_MAX=data.frame(resMAX)

# Inspect the result tables of significantly upregulated and downregulated genes, at a 5% false discovery rate (FDR) as follows:
resSigMAX=resMAX[which(resMAX$padj<0.05),]
head(resSigMAX[order(resSigMAX$log2FoldChange, decreasing = TRUE), ])
head(resSigMAX[order(resSigMAX$log2FoldChange, decreasing = FALSE), ])

write.csv(resSigMAX[order(resSigMAX$log2FoldChange, decreasing = TRUE), ],file = "DE_max_all.csv")
deg_max=rownames(resSigMAX)
# Count the number of genes with significant differential expression at a FDR of 5%:
table(resMAX$padj < 0.05)

#plot heatmap of DEGs
library(pheatmap)
select_max <- deg_max
nt_max <- normTransform(ddsMAX, f = log2, pc = 1) 
log2.norm.counts_max <- assay(nt_max)[select_max,]
df_max <- colData_max
sampleinfo.max=data.frame(group=samples_max$group)
rownames(sampleinfo.max)=colnames(log2.norm.counts_max)
ann_colors_max=list(group=c(Max_T="red", Max_C="grey"))
heatmap_max=pheatmap(log2.norm.counts_max, 
                     cluster_rows=TRUE, 
                     show_rownames=TRUE,
                     cluster_cols=TRUE,
                     scale = "row",
                     border_color = NA,
                     annotation_col=sampleinfo.max,
                     annotation_colors = ann_colors_max,
                     annotation_names_col = F)
heatmap_max


##perform GO enrichment analysis 
#perform GO analysis on MIN
library(topGO)
#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get go id for universe
en_genes=rownames(ddsMAX)
combine=read.csv("combine_within.csv",header = T)
combine=combine[,-1]
en_genes=combine[combine$ga %in% en_genes,]$xm
gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id","go_id"), 
                values = en_genes, 
                mart = mart)

#get all go id for DE max genes
DE_MAX_star=data.frame(deg_max)
en_genes_DE_max=combine[combine$ga %in% DE_MAX_star$deg_max,]
en_genes_DE_max=en_genes_DE_max$xm
DE.gene.max=getBM(filters = "ensembl_gene_id", 
                   attributes = c("ensembl_gene_id","go_id"), 
                   values = en_genes_DE_max, 
                   mart = mart) 

# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
DE.gene.max <- DE.gene.max[DE.gene.max$go_id != '',] 

# convert from table format to list format
geneID2GO.max <- by(gene_pool$go_id,
                     gene_pool$ensembl_gene_id,
                     function(x) as.character(x))

myInterestingGenes.max=by(DE.gene.max$go_id,
                           DE.gene.max$ensembl_gene_id,
                           function(x) as.character(x))

# examine result
head(geneID2GO.max)
head(myInterestingGenes.max)

correction<-"fdr"
geneNames.max = names(geneID2GO.max)
myInterestingGenesNames.max=names(myInterestingGenes.max)
geneList.max = factor(as.integer(geneNames.max %in% myInterestingGenesNames.max))
names(geneList.max) <- geneNames.max

ontology=c("MF","BP","CC")
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList.max, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO.max)
  fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
  fisherResCor = p.adjust(score(fisherRes), method=correction)
  weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
  weightResCor = p.adjust(score(weightRes), method=correction)
  allRes = GenTable(tgData, 
                    classic=fisherRes, 
                    weight=weightRes, 
                    orderBy="weight", 
                    ranksOf="classic", 
                    topNodes=200)
  allRes$fisher.COR = fisherResCor[allRes$GO.ID]
  allRes$weight.COR = weightResCor[allRes$GO.ID]
  write.csv(allRes, paste("DEgenes_max_star",ontology[i],"csv",sep="."))
} 
