library(data.table)
library(dplyr)
library(readr)

args=commandArgs(TRUE)

ann=args[1]
linx=args[2]

sv_ann=fread(ann)
colnames(sv_ann)=c("CHROM","POS","ID","QUAL","LENGTH","TYPE","PURPLE_VAF","PURPLE_CN")
sv_linx=fread(linx)

sv_ann$PURPLE_VAF=sapply(strsplit(sv_ann$PURPLE_VAF, ","), function(x) x[1])
sv_ann$PURPLE_CN=sapply(strsplit(sv_ann$PURPLE_CN, ","), function(x) x[1])
sv_ann[]=lapply(sv_ann, function(x) replace(x, x == ".", NA_character_))

sv_ann$SvID=NA_character_
sv_ann$ClusterID=NA_character_
sv_ann$ClusterReason=NA_character_
sv_ann$Clustered=NA_character_
sv_ann$Foldback=NA_character_
sv_ann$FragileSite=NA_character_
sv_ann$LineType=NA_character_
sv_ann$LocalTopology=NA_character_

for (i in 1:nrow(sv_ann)){
    if (sv_ann$ID[i] %in% sv_linx$vcfId){
      sv_ann$SvID[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"svId"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"svId"])),NA_character_)
      sv_ann$ClusterID[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"clusterId"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"clusterId"])),NA_character_)
      sv_ann$ClusterReason[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"clusterReason"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"clusterReason"])),NA_character_)
      sv_ann$Foldback[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"isFoldback"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"isFoldback"])),NA_character_)
      sv_ann$FragileSite[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"fragileSiteStart"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"fragileSiteStart"])),NA_character_)
      sv_ann$LineType[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"lineTypeStart"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"lineTypeStart"])),NA_character_)
      sv_ann$LocalTopology[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"localTopologyStart"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==sv_ann$ID[i] ),"localTopologyStart"])),NA_character_)
    }else if(paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o") %in% sv_linx$vcfId){
      sv_ann$SvID[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"svId"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"svId"])),NA_character_)
      sv_ann$ClusterID[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"clusterId"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"clusterId"])),NA_character_)
      sv_ann$ClusterReason[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"clusterReason"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o") ),"clusterReason"])),NA_character_)
      sv_ann$Foldback[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"isFoldback"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"isFoldback"])),NA_character_)
      sv_ann$FragileSite[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"fragileSiteEnd"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"fragileSiteEnd"])),NA_character_)
      sv_ann$LineType[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"lineTypeEnd"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"lineTypeEnd"])),NA_character_)
      sv_ann$LocalTopology[i]=ifelse(nzchar(unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"localTopologyEnd"]))),unname(unlist(sv_linx[which(sv_linx$vcfId==paste0(substring(sv_ann$ID[i], 1, nchar(sv_ann$ID[i])-1), "o")),"localTopologyEnd"])),NA_character_)
    }
  sv_ann$Clustered[i]=ifelse(is.na(sv_ann$ClusterReason[i]),FALSE,TRUE)
  }
  
readr::write_tsv(sv_ann,file = strsplit(ann, "/")[[1]][length(strsplit(ann, "/")[[1]])])
