###########################################################
## Differentially Methylated Sites use MAX as an example ##
###########################################################

library(methylKit)
library(graphics)
library(GenomicRanges)
library(IRanges)
library(mixtools)
library(data.table)


filenames=list.files(path=directory, full.names=TRUE)
filenames=as.list(filenames)
names=list.files(path=directory)
names=gsub(".cov.gz", "", names)
names=as.list(names)

# Metadata for all samples
info=read.csv("samples_ren.csv")
library(dplyr)
info=info %>%
  slice(match(names, libraryname))

my.methRaw=methRead(location = filenames,
                    sample.id = names,
                    assembly = "9spine",
                    pipeline = 'bismarkCoverage',
                    context = "CpG",
                    treatment = c(0,0,0,1,1,1,0,0,0,1,1,1), # control is 0, treatment is 1
                    mincov = 5)

# filter dataset for minimum coverage = 5 and remove the 99.9th percentile
filtered.my.methRaw = filterByCoverage(my.methRaw, 
                                       lo.count=5,
                                       lo.perc = NULL,
                                       hi.count = NULL,
                                       hi.perc = 99.9)

# normalize read coverages between samples to avoid bias introduced by systematically more sequenced samples

normalized.myobj=normalizeCoverage(filtered.my.methRaw, method="median")

# merging samples (step should be saved)------

meth.all=unite(normalized.myobj, destrand = F, mc.cores = 8)

save(meth.all, info, file = "diffmeth_mosquitofish.RData")

# load("./diffmeth_mosquitofish.RData")
##########################################
#Exclude SNPs for methylation analyses####
##########################################
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

# The C/T and A/G location is produced before purging high LD sites, using parental samples
# load CT SNPs and duplicate column and save as bedfile, result_pat_C_T_maf.bed
CT_maf= read.csv(file="./result_pat_C_T.txt", sep="\t", header=FALSE)
CT_maf = cbind(CT_maf,V3=rep(CT_maf$V2))
write.table(CT_maf,file="./result_pat_C_T_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_CT_bed<- bed_to_granges(file="./result_pat_C_T_mafSTARTEND.bed")

# load GA SNPs and duplicate column and save as bedfile, result_pat_G_A_maf.bed is created in step 5.2
AG_maf= read.csv(file="./result_pat_G_A.txt", sep="\t", header=FALSE)
AG_maf = cbind(AG_maf,V3=rep(AG_maf$V2))
write.table(AG_maf,file="./result_pat_G_A_mafSTARTEND.bed",row.names = FALSE,quote = FALSE,sep="\t", col.names = FALSE)

blacklist_GA_bed<- bed_to_granges(file="./result_pat_G_A_mafSTARTEND.bed")

# interesect bedfile and unite-file
#### create overlap --> these positions are corrected for CT SNPs
unite_norm_5x_GRanges <- as(meth.all, "GRanges")
Overlap_CT=unite_norm_5x_GRanges[countOverlaps(unite_norm_5x_GRanges, blacklist_CT_bed) == 0L]
meth.1=makeMethylDB(meth.all,"methylBaseDB")
unite_norm_5x_CT <- selectByOverlap(meth.1, granges(Overlap_CT))

# make methylkit database
objDB_CT=makeMethylDB(unite_norm_5x_CT,"methylBaseDB")

##### create overlap --> these positions are corrected for GA SNPs
unite_norm_5x_CT_GRanges <- as(unite_norm_5x_CT,"GRanges")
Overlap_GA=unite_norm_5x_CT_GRanges[countOverlaps(unite_norm_5x_CT_GRanges, blacklist_GA_bed) == 0L]
unite_norm_5x_CT_GA <- selectByOverlap(objDB_CT, Overlap_GA)



# Extract max samples for max differential methylation analysis
meth.MAX=reorganize(unite_norm_5x_CT_GA, sample.ids = c("maxc1","maxc2","maxc3","maxt1","maxt2","maxt3"), treatment = c(0,0,0,1,1,1))

# Calculate differential methylation between treatment and control for each CpG site 
myDiff.MAX=calculateDiffMeth(meth.MAX, mc.cores = 2)

# select differentially methylated bases based on q-value-----
## get hyper methylated bases
myDiff0p.hyper.MAX=getMethylDiff(myDiff.MAX,difference=0,qvalue=0.05,type="hyper")

# get hypo methylated bases
myDiff0p.hypo.MAX=getMethylDiff(myDiff.MAX,difference=0,qvalue=0.05,type="hypo")

# get all differentially methylated bases, and the number of differential methlylated bases
myDiff0p.MAX=getMethylDiff(myDiff.MAX,difference=0,qvalue=0.05)
nrow(myDiff0p.MAX) 
write.table(myDiff0p.MAX,file = "diffmeth_max_diff0.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# For heatmap ploting purpose
# get methylation percentage for all CpG sites
DMCs.max=regionCounts(meth.MAX, regions = as(myDiff0p.MAX, "GRanges"))
perc.DMCs.max=percMethylation(DMCs.max)

# heatmap based on DMCs
library(pheatmap)
sample.max=info[grep("max",info$libraryname),]
sampleinfo.max=data.frame(group=sample.max$group)
rownames(sampleinfo.max)=colnames(perc.DMCs.max)
ann_colors=list(group=c(Max_T="red", Max_C="grey"))
p.heatmap.max=pheatmap(perc.DMCs.max, 
                       cluster_rows=TRUE, 
                       show_rownames=FALSE,
                       show_colnames = TRUE,
                       cluster_cols=TRUE,
                       border_color = NA,
                       scale = "row",
                       clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean",
                       clustering_method = "complete",
                       annotation_col = sampleinfo.max,
                       annotation_colors = ann_colors,
                       annotation_names_col = F)
p.heatmap.max



# Identify gene associated with DMCs-------
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPpeakAnno)

mosquitofish=makeTxDbFromGFF("./mosquitofish.gff3",
                             format = "gff3")
mosquitofish.gene=genes(mosquitofish)

scaffold=seqlevels(mosquitofish.gene)
myDiff0p.MAX.grange=as(myDiff0p.MAX,"GRanges")

###rename seqlevels to satisfy the need in annotatePeakInBatch()
for (i in 1:length(scaffold)) {
  seqlevels(myDiff0p.MAX.grange)=gsub(scaffold[i],i,seqlevels(myDiff0p.MAX.grange))
}

for (i in 1:length(scaffold)) {
  seqlevels(mosquitofish.gene)=gsub(scaffold[i],i,seqlevels(mosquitofish.gene))
}

overlaps.dmc.gene.max=annotatePeakInBatch(myDiff0p.MAX.grange, 
                                          AnnotationData = mosquitofish.gene,
                                          output = "overlapping",
                                          maxgap = 0)

a=data.frame(overlaps.dmc.gene.max)

# remove no feature for some of the scaffold. 
a=a[!is.na(a$feature),]

length(unique(a$feature)) 
dmc_max=a[!duplicated(a$feature),]
write.csv(dmc_max,"dmc_genes_max_0p.csv")

#GO enrichment analysis
#perform GO analysis on MAX
library(topGO)
#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

combine=read.csv("combine_within.csv",sep = ",")
combine=combine[,-1]

dmc_genes_max=read.csv("dmc_genes_max_0p.csv",header = T)
dmc_genes_max=dmc_genes_max$feature


#get go ID for universe
unite_norm_5x_CT_GA_grange=as(unite_norm_5x_CT_GA,"GRanges")

for (i in 1:length(scaffold)) {
  seqlevels(unite_norm_5x_CT_GA_grange)=gsub(scaffold[i],i,seqlevels(unite_norm_5x_CT_GA_grange))
}

pool_gene=annotatePeakInBatch(unite_norm_5x_CT_GA_grange, 
                              AnnotationData = mosquitofish.gene,
                              output = "overlapping",
                              maxgap = 0) 
pool_gene=data.frame(pool_gene)
pool_gene=pool_gene[!is.na(pool_gene$fromOverlappingOrNearest),]
pool_gene=unique(pool_gene$feature)
en_gene=combine[combine$ga %in% pool_gene,]$xm
gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id","go_id"), 
                values = en_gene, 
                mart = mart)

#get all go id for dmc associated genes
en_genes_dmc_max=combine[combine$ga %in% dmc_genes_max,]$xm

dmc.genes.max=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = en_genes_dmc_max, 
                    mart = mart)


# Remove blank entries
gene_pool <- gene_pool[gene_pool$go_id != '',]
dmc.genes.max <- dmc.genes.max[dmc.genes.max$go_id != '',] 

# convert from table format to list format
geneID2GO.max <- by(gene_pool$go_id,
                    gene_pool$ensembl_gene_id,
                    function(x) as.character(x))

myInterestingGenes.max=by(dmc.genes.max$go_id,
                          dmc.genes.max$ensembl_gene_id,
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
  write.csv(allRes, paste("diff0_dmc_genes_max",ontology[i],"csv",sep="."))
} 
