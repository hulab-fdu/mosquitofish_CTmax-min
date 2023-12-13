######################## PCA for alternative splicing analysis using MAX as an example##################################

library(reshape)

#header of five files have been removed and renamed as "_all_fixed.txt"

#A3SS
A3SS=read.table("A3SS_all_fixed.txt",header=FALSE)
colnames(A3SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A3SS=transform(A3SS,IJC=colsplit(IJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
A3SS=transform(A3SS,SJC=colsplit(SJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
A3SS=transform(A3SS,IncLevel=colsplit(IncLevel,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))

A3SS$IJC_sum=rowSums(A3SS$IJC)
A3SS$SJC_sum=rowSums(A3SS$SJC)
A3SS$IncLevel_mean=rowMeans(A3SS$IncLevel,na.rm = TRUE)
write.csv(A3SS,file="A3SS_all.csv")
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
write.csv(A3SS_count_over_20,file="A3SS_all_count_over_20.csv")

#A5SS
A5SS=read.table("A5SS_all_fixed.txt",header=FALSE)
colnames(A5SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A5SS=transform(A5SS,IJC=colsplit(IJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
A5SS=transform(A5SS,SJC=colsplit(SJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
A5SS=transform(A5SS,IncLevel=colsplit(IncLevel,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))

A5SS$IJC_sum=rowSums(A5SS$IJC)
A5SS$SJC_sum=rowSums(A5SS$SJC)
A5SS$IncLevel_mean=rowMeans(A5SS$IncLevel,na.rm = TRUE)
write.csv(A5SS,file="A5SS_all.csv")
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
write.csv(A5SS_count_over_20,file="A5SS_all_count_over_20.csv")


#SE
SE=read.table("SE_all_fixed.txt",header=FALSE)
colnames(SE)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
SE=transform(SE,IJC=colsplit(IJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
SE=transform(SE,SJC=colsplit(SJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
SE=transform(SE,IncLevel=colsplit(IncLevel,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))

SE$IJC_sum=rowSums(SE$IJC)
SE$SJC_sum=rowSums(SE$SJC)
SE$IncLevel_mean=rowMeans(SE$IncLevel,na.rm = TRUE)
write.csv(SE,file="SE_all.csv")
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,] 
write.csv(SE_count_over_20,file="SE_all_count_over_20.csv")


#MXE
MXE=read.table("MXE_all_fixed.txt",header=FALSE)
colnames(MXE)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2ndExonStart_0base","2ndExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
MXE=transform(MXE,IJC=colsplit(IJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
MXE=transform(MXE,SJC=colsplit(SJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
MXE=transform(MXE,IncLevel=colsplit(IncLevel,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))

MXE$IJC_sum=rowSums(MXE$IJC)
MXE$SJC_sum=rowSums(MXE$SJC)
MXE$IncLevel_mean=rowMeans(MXE$IncLevel,na.rm = TRUE)
write.csv(MXE,file="MXE_all.csv")
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
write.csv(MXE_count_over_20,file="MXE_all_count_over_20.csv") 

#RI
RI=read.table("RI_all_fixed.txt",header=FALSE)
colnames(RI)=c("ID","GeneID","geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
RI=transform(RI,IJC=colsplit(IJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
RI=transform(RI,SJC=colsplit(SJC,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))
RI=transform(RI,IncLevel=colsplit(IncLevel,split=",",names=c("maxt1","maxt2","maxt3","maxc1","maxc2","maxc3")))

RI$IJC_sum=rowSums(RI$IJC)
RI$SJC_sum=rowSums(RI$SJC)
RI$IncLevel_mean=rowMeans(RI$IncLevel,na.rm = TRUE)
write.csv(RI,file="RI_all.csv")
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
write.csv(RI_count_over_20,file="RI_all_count_over_20.csv") 


# Create colomn data to plot all filtered splicing events across all clean samples.
A3SS_count_over20=read.csv(file="A3SS_all_count_over_20.csv")
A5SS_count_over20=read.csv(file="A5SS_all_count_over_20.csv")
MXE_count_over20=read.csv(file="MXE_all_count_over_20.csv")
RI_count_over20=read.csv(file="RI_all_count_over_20.csv")
SE_count_over20=read.csv(file="SE_all_count_over_20.csv")
A3SS_count_over20$GeneID <- paste("A3SS", substring(A3SS_count_over20$GeneID,6), sep="_") 
A5SS_count_over20$GeneID <- paste("A5SS", substring(A5SS_count_over20$GeneID,6), sep="_")
MXE_count_over20$GeneID <- paste("MXE", substring(MXE_count_over20$GeneID,6), sep="_")
RI_count_over20$GeneID <- paste("RI", substring(RI_count_over20$GeneID,6), sep="_")
SE_count_over20$GeneID <- paste("SE", substring(SE_count_over20$GeneID,6), sep="_")
all_count_over_20=rbind(A3SS_count_over20[,30:35],A5SS_count_over20[,30:35],MXE_count_over20[,32:37],SE_count_over20[,30:35],RI_count_over20[,30:35]) #substring all Inclevel of each sample and rbind all types of events
colnames(all_count_over_20)=substring(colnames(all_count_over_20),10) #remove "Inclevel."
all_count_over_20_no_na=all_count_over_20[complete.cases(all_count_over_20),] 



# Plot PCA
max_splicing=prcomp(t(all_count_over_20_no_na),center = T)
summary(max_splicing) 
max_splicing_df=as.data.frame(max_splicing$x)
samples=read.csv("samples_ren.csv")
group=as.factor(samples[samples$condition=="MAX",]$group)

#visulize PCA results
g1.max=ggbiplot(max_splicing, 
                obs.scale = 1, 
                var.scale = 1,
                varname.size = 0,
                var.axes = FALSE)
g2.max=g1.max+geom_point(aes(colour = group, fill=group), size=1.5)+
  scale_color_manual(values = c("grey", "red"))
g3.max=g2.max + theme_bw()
print(g3.max)

