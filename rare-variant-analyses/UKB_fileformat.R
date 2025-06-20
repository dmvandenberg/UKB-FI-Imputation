#Script to produce the setlist and annotation file for gene-based testing in Regenie

#example for PTV variants
for (i in 1:22){
  print(i)
  file.name1 <- paste("HC-PTV-chr",i,".",i,".STAAR.variants_table.tsv",sep="")
  file.name2 <- paste("chr",i,"_setlist.txt",sep="")
  file.name3 <- paste("chr",i,"_annotation.txt",sep="")
  f <- read.table(file.name1,header=T,stringsAsFactors = F)
  genes <- unique(f[,4])
  fdat <- NULL
  for (j in 1:length(genes)){
    print(j)
    indg <- which(f[,4]==genes[j])
    vars <- paste(f[indg,1], collapse=",")
    dat <- c(f[indg,4][1], f[indg,2][1], f[indg,3][1], vars)
    fdat <- rbind(fdat, dat)}
  write.table(fdat,file.name2,quote=F, row.names=F, col.names=F, sep="\t")
  fann <- f[,c(1,4)]
  fann$type <- rep("PTV",length(fann[,1]))
  write.table(fann,file.name3,quote=F, row.names=F, col.names=F, sep="\t")
}

#example for missense variants and filter missense variant based on MPC score only 

for (i in 22:22){
  print(i)
  file.name1 <- paste("MIS-chr",i,".tsv",sep="")
  file.name2 <- paste("chr",i,"_missetlist.txt",sep="")
  file.name3 <- paste("chr",i,"_misannotations.txt",sep="")
  f <- read.table(file.name1,header=T,stringsAsFactors = F)
  genes <- unique(f[,4])
  fdat <- NULL
  for (j in 1:length(genes)){
    print(j)
    indg <- which(f[,4]==genes[j])
    vars <- paste(f[indg,1], collapse=",")
    dat <- c(f[indg,4][1], f[indg,2][1], f[indg,3][1], vars)
    fdat <- rbind(fdat, dat)}
  write.table(fdat,file.name2,quote=F, row.names=F, col.names=F, sep="\t")
  fann1 <- f[f$score<2,]
  fann2 <- f[f$score<3&f$score>=2,]
  fann3 <- f[f$score>3,]
  fann1 <- fann1[,c(1,4)]
  fann2 <- fann2[,c(1,4)]
  fann3 <- fann3[,c(1,4)]
  fann1$type <- rep("MPC1",length(fann1[,1]))
  fann2$type <- rep("MPC2",length(fann2[,1]))
  fann3$type <- rep("MPC3",length(fann3[,1]))
  fann <- rbind(fann1,fann2,fann3)
  write.table(fann,file.name3,quote=F, row.names=F, col.names=F, sep="\t")
}

#example to calculate the number of variants for each gene. 
#the results can be used to calculate the number of rare varinats in each individual for downstream analysis.

library(Matrix)
for (i in 1:22){
  print(i)
  file.name1 <- paste("HC-PTV-chr",i,".",i,".STAAR.samples_table.tsv",sep="")
  file.name2 <- paste("HC-PTV-chr",i,".",i,".STAAR.variants_table.tsv",sep="")
  file.name3 <- paste("chr",i,"_ptv_gnomad_af0.1.txt",sep="")
  matrix <- readRDS(file.name1)
  summ <- summary(matrix)
  count <- data.frame(eid=rownames(matrix)[summ$i],variant=colnames(matrix)[summ$j],count=summ$x)
  info <- read.table(file.name2,header=T)
  info <- info[,c("varID","chrom","ENST")]
  colnames(info)[1] <- "variant"
  count <- merge(count,info,by="variant")
  count <- count[,c("eid","count","chrom","ENST")]
  write.table(count,file.name3,row.names=F)
}