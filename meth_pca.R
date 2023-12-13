######PCA for methylation pattern using MAX as an example######

library(methylKit)

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

# normalize read coverages between samples to avoid bias introduced by systematically more sequnenced sameples

normalized.myobj=normalizeCoverage(filtered.my.methRaw, method="median")

# merging samples (step should be saved)------

meth.all=unite(normalized.myobj, destrand = F, mc.cores = 8)

save(meth.all, info, file = "meth_mosquitofish.RData")

# load("./meth.RData")
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
Overlap_GA=unite_norm_5x_CT_GRanges[countOverlaps(unite_norm_5x_CT_GRanges, blacklist_GA_bed) == 0L] #strand
unite_norm_5x_CT_GA <- selectByOverlap(objDB_CT, Overlap_GA)


# PCA----
library(ggbiplot)
library(factoextra)

meth.count=regionCounts(unite_norm_5x_CT_GA, regions = as(unite_norm_5x_CT_GA, "GRanges")) 
perc.meth.all=percMethylation(meth.count)

perc.meth.max=perc.meth.all[,colnames(perc.meth.all)==info[info$condition=="MAX",]$libraryname]

# check if the order of meth count table is the same as info
colnames(perc.meth.max)==info[info$condition=="MAX",]$libraryname

pca.max=prcomp(t(perc.meth.max), center = T)
summary(pca.max)

#Define factors
info.max=info[info$condition=="MAX",]
group=as.factor(info.max[info.max$condition=="MAX",]$group)

#visulize PCA results
g1.max=ggbiplot(pca.max, 
                obs.scale = 1, 
                var.scale = 1,
                varname.size = 0,
                var.axes = FALSE)
g2.max=g1.max+geom_point(aes(colour = group, fill=group), size=4)+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = c("grey","red"))
g3.max=g2.max + theme_bw()
print(g3.max)
