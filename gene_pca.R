library(DESeq2)
library(ggbiplot)
setwd("~/Desktop/htseq_count_star")

#Read in sample table and trim count table
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

keep_gene=rowSums(counts(dds)>=1)>=nrow(colData) # remove genes with less than 1 reads per sample in all samples, note in Huang et al. paper, they only filter gene in half of samples
dds=dds[keep_gene,]
gene_rlog=rlog(dds, blind = F)
gene_normalized=data.frame(gene_rlog@assays@data@listData)

#Trim sample table for CTmax and control
keep_samples_max=grep("MAX", samples_star$condition)
samples_max=samples_star[keep_samples_max,]

#Trim counts table for DE analysis
max_keep=grep("max", colnames(counts_star)) 
max=gene_normalized[,max_keep]

#perform PCA on max samples
variable=samples_max$group

colnames(max)==as.character(samples_max$libraryname) #check the order

#perform PCA 
pca_max=prcomp(t(max), center = T)

#define factors
group=as.factor(variable)

#visualize PCA results
g1=ggbiplot(pca_max, 
            obs.scale = 1, 
            var.scale = 1,
            varname.size=0,
            var.axes = FALSE,
)
g2=g1+geom_point(aes(colour=group, fill=group), size=3)+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("grey", "red"))+
  coord_equal(ratio = 1.5)
g3=g2 + theme_bw()
print(g3)

