######################## DSG analysis for CTmax as an example##################################

library("reshape")
setwd("~/Desktop/rmats_check/")

#header of five files have been removed and renamed as "_all_fixed.txt"

#A3SS
A3SS_max=read.table("A3SS_all_fixed.txt",header=FALSE)
colnames(A3SS_max)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
A3SS_max=transform(A3SS_max,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
A3SS_max=transform(A3SS_max,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
A3SS_max=transform(A3SS_max,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
A3SS_max=transform(A3SS_max,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
A3SS_max=transform(A3SS_max,IncLevel1=colsplit(IncLevel1,split=",",names=c('maxt1','maxt2','maxt3')))
A3SS_max=transform(A3SS_max,IncLevel2=colsplit(IncLevel2,split=",",names=c('maxc1','maxc2','maxc3')))

nrow(A3SS_max)
A3SS_max$IJC_sum=rowSums(A3SS_max$IJC_SAMPLE_1)+rowSums(A3SS_max$IJC_SAMPLE_2)
A3SS_max$SJC_sum=rowSums(A3SS_max$SJC_SAMPLE_1)+rowSums(A3SS_max$SJC_SAMPLE_2)
A3SS_max_count_over_20=A3SS_max[A3SS_max$IJC_sum>20&A3SS_max$SJC_sum>20,] #confident alternative splicing events

A3SS_max_sig=A3SS_max_count_over_20[A3SS_max_count_over_20$FDR<0.05,] #DS events
nrow(A3SS_max_sig) 

#A5SS
A5SS_max=read.table("A5SS_all_fixed.txt",header=FALSE)
colnames(A5SS_max)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
A5SS_max=transform(A5SS_max,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
A5SS_max=transform(A5SS_max,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
A5SS_max=transform(A5SS_max,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
A5SS_max=transform(A5SS_max,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
A5SS_max=transform(A5SS_max,IncLevel1=colsplit(IncLevel1,split=",",names=c('maxt1','maxt2','maxt3')))
A5SS_max=transform(A5SS_max,IncLevel2=colsplit(IncLevel2,split=",",names=c('maxc1','maxc2','maxc3')))

nrow(A5SS_max)
A5SS_max$IJC_sum=rowSums(A5SS_max$IJC_SAMPLE_1)+rowSums(A5SS_max$IJC_SAMPLE_2)
A5SS_max$SJC_sum=rowSums(A5SS_max$SJC_SAMPLE_1)+rowSums(A5SS_max$SJC_SAMPLE_2)
A5SS_max_count_over_20=A5SS_max[A5SS_max$IJC_sum>20&A5SS_max$SJC_sum>20,] 

A5SS_max_sig=A5SS_max_count_over_20[A5SS_max_count_over_20$FDR<0.05,]
nrow(A5SS_max_sig) 


#SE
SE_max=read.table("SE_all_fixed.txt",header=FALSE)
colnames(SE_max)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
SE_max=transform(SE_max,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
SE_max=transform(SE_max,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
SE_max=transform(SE_max,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
SE_max=transform(SE_max,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
SE_max=transform(SE_max,IncLevel1=colsplit(IncLevel1,split=",",names=c('maxt1','maxt2','maxt3')))
SE_max=transform(SE_max,IncLevel2=colsplit(IncLevel2,split=",",names=c('maxc1','maxc2','maxc3')))


nrow(SE_max)
SE_max$IJC_sum=rowSums(SE_max$IJC_SAMPLE_1)+rowSums(SE_max$IJC_SAMPLE_2)
SE_max$SJC_sum=rowSums(SE_max$SJC_SAMPLE_1)+rowSums(SE_max$SJC_SAMPLE_2)
SE_max_count_over_20=SE_max[SE_max$IJC_sum>20&SE_max$SJC_sum>20,]

SE_max_sig=SE_max_count_over_20[SE_max_count_over_20$FDR<0.05,]
nrow(SE_max_sig) 

#MXE
MXE_max=read.table("MXE_all_fixed.txt",header=FALSE)
colnames(MXE_max)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2ndExonStart_0base","2ndExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
MXE_max=transform(MXE_max,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
MXE_max=transform(MXE_max,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
MXE_max=transform(MXE_max,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
MXE_max=transform(MXE_max,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
MXE_max=transform(MXE_max,IncLevel1=colsplit(IncLevel1,split=",",names=c('maxt1','maxt2','maxt3')))
MXE_max=transform(MXE_max,IncLevel2=colsplit(IncLevel2,split=",",names=c('maxc1','maxc2','maxc3')))

nrow(MXE_max)
MXE_max$IJC_sum=rowSums(MXE_max$IJC_SAMPLE_1)+rowSums(MXE_max$IJC_SAMPLE_2)
MXE_max$SJC_sum=rowSums(MXE_max$SJC_SAMPLE_1)+rowSums(MXE_max$SJC_SAMPLE_2)
MXE_max_count_over_20=MXE_max[MXE_max$IJC_sum>20&MXE_max$SJC_sum>20,] 

MXE_max_sig=MXE_max_count_over_20[MXE_max_count_over_20$FDR<0.05,]
nrow(MXE_max_sig) 

#RI
RI_max=read.table("RI_all_fixed.txt",header=FALSE)
colnames(RI_max)=c("ID","GeneID","geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","lncFormLen","SkipFormLen","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
RI_max=transform(RI_max,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
RI_max=transform(RI_max,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
RI_max=transform(RI_max,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('maxt1','maxt2','maxt3')))
RI_max=transform(RI_max,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('maxc1','maxc2','maxc3')))
RI_max=transform(RI_max,IncLevel1=colsplit(IncLevel1,split=",",names=c('maxt1','maxt2','maxt3')))
RI_max=transform(RI_max,IncLevel2=colsplit(IncLevel2,split=",",names=c('maxc1','maxc2','maxc3')))

nrow(RI_max)
RI_max$IJC_sum=rowSums(RI_max$IJC_SAMPLE_1)+rowSums(RI_max$IJC_SAMPLE_2)
RI_max$SJC_sum=rowSums(RI_max$SJC_SAMPLE_1)+rowSums(RI_max$SJC_SAMPLE_2)
RI_max_count_over_20=RI_max[RI_max$IJC_sum>20&RI_max$SJC_sum>20,] 

RI_max_sig=RI_max_count_over_20[RI_max_count_over_20$FDR<0.05,]
nrow(RI_max_sig)


## AS events
A3SS_max_count_over_20$Type=rep("A3SS",nrow(A3SS_max_count_over_20))
A5SS_max_count_over_20$Type=rep("A5SS",nrow(A5SS_max_count_over_20))
MXE_max_count_over_20$Type=rep("MXE",nrow(MXE_max_count_over_20))
SE_max_count_over_20$Type=rep("SE",nrow(SE_max_count_over_20))
RI_max_count_over_20$Type=rep("RI",nrow(RI_max_count_over_20))

AS_max=rbind(rbind(rbind(rbind(A3SS_max_count_over_20[,c(2,20,23,26)]
                               ,A5SS_max_count_over_20[,c(2,20,23,26)])
                         ,MXE_max_count_over_20[,c(2,22,25,28)])
                   ,SE_max_count_over_20[,c(2,20,23,26)])
             ,RI_max_count_over_20[,c(2,20,23,26)])
AS_max$absIncleveldifference=abs(AS_max$IncLevelDifference)
write.csv(AS_max,file="AS_events_max.csv") # Total confident AS events between treatment and control
AS_genes_max=AS_max[!duplicated(AS_max$GeneID),] 

####DSG
A3SS_max_sig$Type=rep("A3SS",nrow(A3SS_max_sig))
A5SS_max_sig$Type=rep("A5SS",nrow(A5SS_max_sig))
MXE_max_sig$Type=rep("MXE",nrow(MXE_max_sig))
SE_max_sig$Type=rep("SE",nrow(SE_max_sig))
RI_max_sig$Type=rep("RI",nrow(RI_max_sig))

DSG_max=rbind(rbind(rbind(rbind(A3SS_max_sig[,c(2,20,23,26)]
                                ,A5SS_max_sig[,c(2,20,23,26)])
                          ,MXE_max_sig[,c(2,22,25,28)])
                    ,SE_max_sig[,c(2,20,23,26)])
              ,RI_max_sig[,c(2,20,23,26)])
DSG_max$absIncleveldifference=abs(DSG_max$IncLevelDifference)

# The DS event with the largest absIncleveldifference of a DSG represent the differential splicing level of the gene.
DSG_max_unique_gene=DSG_max[order(DSG_max[,'GeneID'],-DSG_max[,'absIncleveldifference']),]
DSG_max_unique_gene=DSG_max_unique_gene[!duplicated(DSG_max_unique_gene$GeneID),]
write.csv(DSG_max_unique_gene,file="DSG_max_unique.csv") # Total unique DS genes between treatment and control


#GO enrichment analysis
combine=read.csv("combine.csv",header = TRUE)
combine=combine[,-1]

library(topGO)
#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))
combine=read.csv("combine.csv",header=T)
combine=combine[,-1]

#get all go id for universe (AS_genes)
en_genes=combine[combine$ga %in% AS_genes_max$GeneID,]$xm
gene_pool=getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id","go_id"), 
                values = en_genes, 
                mart = mart)

#get all go id for DS genes
DS_genes_max=data.frame(DSG_max_unique_gene$GeneID)
en_genes_DS_max=combine[combine$ga %in% DS_genes_max$DSG_max_unique_gene.GeneID,]
en_genes_DS_max=en_genes_DS_max$xm
DS.gene.max=getBM(filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id","go_id"), 
                  values = en_genes_DS_max, 
                  mart = mart) 

# Remove blank entries
gene_pool<- gene_pool[gene_pool$go_id != '',]
DS.gene.max <- DS.gene.max[DS.gene.max$go_id != '',] 

# convert from table format to list format
geneID2GO <- by(gene_pool$go_id,
                gene_pool$ensembl_gene_id,
                function(x) as.character(x))

myInterestingGenes=by(DS.gene.max$go_id,
                      DS.gene.max$ensembl_gene_id,
                      function(x) as.character(x))

# examine result
head(geneID2GO)
head(myInterestingGenes)

correction<-"fdr"
geneNames = names(geneID2GO)
myInterestingGenesNames=names(myInterestingGenes)
geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))
names(geneList) <- geneNames

ontology=c("MF","BP","CC")
for (i in 1:length(ontology)) {
  tgData = new("topGOdata", 
               ontology = ontology[i], 
               allGenes = geneList, 
               annot = annFUN.gene2GO, 
               gene2GO = geneID2GO)
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
  write.csv(allRes, paste("DSgenes_max",ontology[i],"csv",sep="."))
} 

