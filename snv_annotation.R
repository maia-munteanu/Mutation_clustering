library(data.table)
library(VariantAnnotation)
library(dplyr)
library(stringi)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)
library(ggplot2)
library(gtools)
library(pals)
library(ggnewscale)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
dnaRevCompl <- function(nucSeq){return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))}

args=commandArgs(TRUE)

cores=as.integer(args[1])
closer=as.integer(args[2])
close=as.integer(args[3])
input_file=fread(args[4])
assembly=args[5]
chr_sizes=fread(args[6]); colnames(chr_sizes)=c("CHROM","SIZE")
name=args[7]
tsv=fread(args[8])
snvsnv=fread(args[9])
sv=fread(args[10])
ratio=as.numeric(args[11])
ocloser=as.numeric(args[12])
oclose=as.numeric(args[13])
ounclustered=as.numeric(args[14])
sizecloser=as.numeric(args[15])
sizeclose=as.numeric(args[16])
sizeunclustered=as.numeric(args[17])
closer_decomp=fread(args[18])
closer_denovo=fread(args[19])
close_decomp=fread(args[20])
close_denovo=fread(args[21])
unclustered_decomp=fread(args[22])
unclustered_denovo=fread(args[23])

#change signature names
colnames(closer_denovo)[3:ncol(closer_denovo)]=paste0(colnames(closer_denovo)[3:ncol(closer_denovo)],"-Closer"); colnames(closer_denovo)[3:ncol(closer_denovo)]=gsub("SBS96","",colnames(closer_denovo)[3:ncol(closer_denovo)])
colnames(close_denovo)[3:ncol(close_denovo)]=paste0(colnames(close_denovo)[3:ncol(close_denovo)],"-Close"); colnames(close_denovo)[3:ncol(close_denovo)]=gsub("SBS96","",colnames(close_denovo)[3:ncol(close_denovo)])
colnames(unclustered_denovo)[3:ncol(unclustered_denovo)]=paste0(colnames(unclustered_denovo)[3:ncol(unclustered_denovo)],"-Unclustered"); colnames(unclustered_denovo)[3:ncol(unclustered_denovo)]=gsub("SBS96","",colnames(unclustered_denovo)[3:ncol(unclustered_denovo)])

#change chrom names to chrX instead of just X (we can switch back at the end)
colnames(tsv)[colnames(tsv) == "#CHROM"] <- "CHROM"
snvsnv=snvsnv[,-c("sample","sbs96")]; colnames(snvsnv)=c("CHROM","POS","REF","ALT","SNV_FDR","SNV_ObsDistance","SNV_ExpDistance","SNV_ClusterType","SNV_ClusterSize","SNV_ClusterID")
if(length(grep("chr",chr_sizes$CHROM))==0){chr_sizes$CHROM=paste0("chr",chr_sizes$CHROM)}
if(length(grep("chr",tsv$CHROM))==0){tsv$CHROM=paste0("chr",tsv$CHROM)}
if(length(grep("chr",snvsnv$CHROM))==0){snvsnv$CHROM=paste0("chr",snvsnv$CHROM)}
if(length(grep("chr",sv$CHROM))==0){sv$CHROM=paste0("chr",sv$CHROM)}

#add absolute position
tsv %>% group_by(CHROM) %>% mutate(POS_ABS=ifelse(CHROM=="chr1", as.numeric(POS), as.numeric(POS)+sum(chr_sizes$SIZE[1:(which(chr_sizes$CHROM == unique(CHROM))-1)]))) %>% pull(POS_ABS) -> tsv$POS_ABS

# add sample info (ratio of randomised to observed, closer/close/uncl counts, closer/close/uncl nt at risk)
tsv$RadomisedSNVRatio=ratio
tsv$ObservedCloser=ocloser
tsv$ObservedClose=oclose
tsv$ObservedUnclustered=ounclustered
tsv$CloserAtRisk=sizecloser
tsv$CloseAtRisk=sizeclose
tsv$UnclusteredAtRisk=sizeunclustered

#get the upstream and downstream base for each variant (assembly used is important here)
if (assembly=="hg19"){tsv$UP=as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19,tsv$CHROM,start=tsv$POS-1, end=tsv$POS-1))$x; tsv$DOWN=as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19,tsv$CHROM,start=tsv$POS+1, end=tsv$POS+1))$x
}else if (assembly=="hg38"){tsv$UP=as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg38,tsv$CHROM,start=tsv$POS-1, end=tsv$POS-1))$x; tsv$DOWN=as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg38,tsv$CHROM,start=tsv$POS+1, end=tsv$POS+1))$x
}else{stop("Assembly provided is not hg19 or hg38")}

#add mutation contexts (64 full or 32 reduced by pyrimidine)
tsv$TNC64=paste0(tsv$UP,tsv$REF,tsv$DOWN)
tsv$TNC32=ifelse(tsv$REF %in% c("C","T"),paste0(tsv$UP,tsv$REF,tsv$DOWN),paste0(dnaRevCompl(tsv$DOWN),dnaRevCompl(tsv$REF),dnaRevCompl(tsv$UP)))

#add mutation (12 full without context, 6 reduced by pyrimidine without context, 192 full with context, 96 reduced by pyrimidine with context)
tsv$MUT12=paste0(tsv$REF,">",tsv$ALT)
tsv$MUT6=ifelse(tsv$REF %in% c("C","T"), paste0(tsv$REF,">",tsv$ALT), paste0(dnaRevCompl(tsv$REF),">",dnaRevCompl(tsv$ALT)))
tsv$MUTwc192=paste0(tsv$UP,"[",tsv$REF,">",tsv$ALT,"]",tsv$DOWN)
tsv$MUTwc96=ifelse(tsv$REF %in% c("C","T"), paste0(tsv$UP,"[",tsv$REF,">",tsv$ALT,"]",tsv$DOWN), paste0(dnaRevCompl(tsv$DOWN),"[",dnaRevCompl(tsv$REF),">",dnaRevCompl(tsv$ALT),"]",dnaRevCompl(tsv$UP)))

#add motif information - APOBEC and AID known preferred motifs
tsv$APOBECmotif=ifelse(tsv$MUTwc96 %in% c("T[C>A]A","T[C>A]T","T[C>G]A","T[C>G]T","T[C>T]A","T[C>T]T"),TRUE,FALSE)
if (assembly=="hg19"){tsv$AIDmotif=ifelse(tsv$REF %in% c("C","G"),ifelse(tsv$REF=="C",as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19,tsv$CHROM,start=tsv$POS-2, end=tsv$POS+1))$x,as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19,tsv$CHROM,start=tsv$POS-1, end=tsv$POS+2))$x),FALSE)
}else if (assembly=="hg38"){tsv$AIDmotif=ifelse(tsv$REF %in% c("C","G"),ifelse(tsv$REF=="C",as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg38,tsv$CHROM,start=tsv$POS-2, end=tsv$POS+1))$x,as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg38,tsv$CHROM,start=tsv$POS-1, end=tsv$POS+2))$x),FALSE)
}else{stop("Assembly provided is not hg19 or hg38")}
tsv[which(tsv$REF=="C" & tsv$AIDmotif %in% c("AACC","AACT","AGCC","AGCT","TACC","TACT","TGCC","TGCT")),"AIDmotif"] <- TRUE
tsv[which(tsv$REF=="G" & tsv$AIDmotif %in% c("AGCA","AGCT","AGTA","AGTT","GGCA","GGCT","GGTA","GGTT")),"AIDmotif"] <- TRUE
tsv[-which(tsv$AIDmotif %in% c(TRUE,FALSE)),"AIDmotif"] <- FALSE; tsv$AIDmotif=as.logical(tsv$AIDmotif)

#add inter-mutation distance to the nearest downstream mutation; if the mutation is the last on a chromosome, count distance to end of chromosome
tsv %>% group_by(CHROM) %>% mutate(IMD = lead(POS, 1) - POS) -> tsv
tsv$IMD[is.na(tsv$IMD)]=unname(unlist(chr_sizes[chr_sizes$CHROM %in% tsv$CHROM[is.na(tsv$IMD)],"SIZE"] - tsv$POS[is.na(tsv$IMD)])) #for the last mutation on a chromosome, IMD is the distance to the end of the chromosome

#add snv-snv cluster information (fdr, observed and expected snv-snv distance, cluster type, size and ID); note cluster size is also used when the mutations don't form a cluster per se
tsv <- left_join(tsv, snvsnv, by = c("CHROM", "POS", "REF", "ALT"))
#chnage cluster type names and also add doublet and multi base substitutions when IMD is 1 (David doesn't account for this)
tsv <- tsv %>% mutate(SNV_ClusterType = ifelse(is.na(SNV_ClusterType), "Non-clustered", SNV_ClusterType), SNV_ClusterType = recode(SNV_ClusterType, "omikli" = "Omikli", "kataegis" = "Kataegis")) 
tsv <- tsv %>% group_by(CHROM) %>% mutate(before = lag(IMD), after = lead(IMD), SNV_ClusterType = case_when(IMD == 1 & (is.na(before) | before != 1) & (is.na(after) | after != 1) ~ "DBS", IMD == 1 & ((before == 1) | (after == 1)) ~ "MBS", TRUE ~ as.character(SNV_ClusterType))) %>% ungroup() %>% select(-c(before,after))

#add sv information (take all neighboring SVs info and merge into one Ann field, similar format to the ANN field)
#also add signatures (de novo and decomposed) as well as their probabilities 
registerDoParallel(cores = cores)
signatures <- foreach(i = 1:nrow(tsv), .combine = 'rbind', .packages = c("dplyr", "tidyverse")) %dopar% {
  if (tsv$SVSNV[i]=="CLOSER"){
    svs=subset(sv, CHROM == tsv$CHROM[i] & abs(POS - tsv$POS[i]) <= closer)
    svs=mutate_all(svs, as.character); svs[is.na(svs)] <- ""
    SVSNV_Ann=paste(apply(svs,1,paste0,collapse = "|"),collapse=",")
    denovo_probs=subset(closer_denovo,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(closer_denovo)]
    decomp_probs=subset(closer_decomp,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(closer_decomp)]
  } else if (tsv$SVSNV[i]=="CLOSE"){
    svs=subset(sv, CHROM == tsv$CHROM[i] & abs(POS - tsv$POS[i]) <= close & abs(POS - tsv$POS[i]) > closer)
    svs=mutate_all(svs, as.character); svs[is.na(svs)] <- ""
    SVSNV_Ann=paste(apply(svs,1,paste0,collapse = "|"),collapse=",")
    denovo_probs=subset(close_denovo,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(close_denovo)]
    decomp_probs=subset(close_decomp,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(close_decomp)]
  } else if (tsv$SVSNV[i]=="UNCLUSTERED"){
    SVSNV_Ann = NA_character_
    denovo_probs=subset(unclustered_denovo,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(unclustered_denovo)]
    decomp_probs=subset(unclustered_decomp,`Sample Names` == name & MutationTypes==tsv$MUTwc96[i])[1,3:ncol(unclustered_decomp)] }
  if(!any(is.na(denovo_probs))){set.seed(i)
    SignatureDenovo=names(denovo_probs)[t(rmultinom(1,1,denovo_probs)) %>% {which(. != 0)}]
    SignatureDenovoProb=round(denovo_probs[[SignatureDenovo]],4)
  } else {SignatureDenovo = NA_character_; SignatureDenovoProb = NA}
  if(!any(is.na(decomp_probs))){set.seed(i)
    SignatureDecomp=names(decomp_probs)[t(rmultinom(1,1,decomp_probs)) %>% {which(. != 0)}]
    SignatureDecompProb=round(decomp_probs[[SignatureDecomp]],4)
  } else {SignatureDecomp = NA_character_; SignatureDecompProb = NA}
  data.frame(SVSNV_Ann = SVSNV_Ann, SignatureDenovo = SignatureDenovo, SignatureDenovoProb = SignatureDenovoProb, SignatureDecomp = SignatureDecomp, SignatureDecompProb = SignatureDecompProb)}
tsv <- cbind(tsv, signatures) # Combine the results with the original tsv data frame

#add sample information such as TMB, microsatellite stability, gender, primary tissue, cancer type etc from the input file (doubles as a metadata )
tsv=cbind(tsv,input_file[which(input_file$sample==name),-c("sample","sv","linx","snv")])

#export tsv
write.table(tsv, file=paste0(name,"_annotated.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

# start configuring plot settings
sv %>% group_by(CHROM) %>% mutate(POS_ABS=ifelse(CHROM=="chr1", as.numeric(POS), as.numeric(POS)+sum(chr_sizes$SIZE[1:(which(chr_sizes$CHROM == unique(CHROM))-1)]))) %>% pull(POS_ABS) -> sv$POS_ABS ##get absolut positions for SVs for plotting
sv$TYPE[which(is.na(sv$TYPE))]="UNK"

for(i in 1:nrow(chr_sizes)){chr_sizes$POS[i]=sum(chr_sizes$SIZE[1:i])}
chr_sizes$POS_LABELS=chr_sizes$POS-(chr_sizes$SIZE/2) #add half chromosome positions for ggplot labels

#set color palettes for SVs, SNVs and signatures
sv_levels=c("INS","DEL","INV","DUP","CTX","UNK")
sv$TYPE <- factor(sv$TYPE, levels = sv_levels)
sv_colors=c('#03BDEF','#E42926','#f27e02','#85cf1f','#762ca3','#828181')   
sv_colors=sv_colors[sv_levels %in% unique(sv$TYPE)]

snv_levels=c("C>A","C>G","C>T","T>A","T>C","T>G")
tsv$MUT6 <- factor(tsv$MUT6, levels = snv_levels)
snv_colors=c("#03BDEF","#454343","#E42926","#CBCACA","#A2CF63","#ECC7C5")
snv_colors=snv_colors[snv_levels %in% unique(tsv$MUT6)]

cluster_levels=c("DBS","MBS","Omikli","Kataegis","Non-clustered")
tsv$SNV_ClusterType <- factor(tsv$SNV_ClusterType, levels = cluster_levels)
cluster_colors=c("#ee6677","#aa3377","#66ccee","#4477aa",'#828181')
cluster_colors=cluster_colors[cluster_levels %in% unique(tsv$SNV_ClusterType)]

svsnv=c("CLOSER","CLOSE","UNCLUSTERED")
tsv$SVSNV <- factor(tsv$SVSNV, levels = svsnv)
svsnv_shapes=c(15,17,19)
svsnv_shapes=svsnv_shapes[svsnv %in% unique(tsv$SVSNV)]
svsnv_colors=c("#851667","#b87ca8","#dfc8d9")
svsnv_colors=svsnv_colors[svsnv %in% unique(tsv$SVSNV)]
legend_shapes=c(rep(svsnv_shapes[1],times=length(unique(tsv$SignatureDenovo[which(tsv$SVSNV=="CLOSER")]))),
                rep(svsnv_shapes[2],times=length(unique(tsv$SignatureDenovo[which(tsv$SVSNV=="CLOSE")]))),
                rep(svsnv_shapes[3],times=length(unique(tsv$SignatureDenovo[which(tsv$SVSNV=="UNCLUSTERED")]))))

denovo_sig_levels=c(colnames(closer_denovo)[3:ncol(closer_denovo)],colnames(close_denovo)[3:ncol(close_denovo)],colnames(unclustered_denovo)[3:ncol(unclustered_denovo)])
decomp_sig_levels=mixedsort(unique(c(colnames(closer_decomp)[3:ncol(closer_decomp)],colnames(close_decomp)[3:ncol(close_decomp)],colnames(unclustered_decomp)[3:ncol(unclustered_decomp)])))
tsv$SignatureDenovo <- factor(tsv$SignatureDenovo, levels = denovo_sig_levels)
tsv$SignatureDecomp <- factor(tsv$SignatureDecomp, levels = decomp_sig_levels)
colors=c(rev(glasbey(n = 32)),rev(kelly(n=22))); colors=colors[(!colors %in% c("#000033","#222222","#F3C300","#A1CAF1","#BE0032","#F1085C","#858567","#00479E"))]
denovo_sig_colors <- colors[1:length(denovo_sig_levels)]
decomp_sig_colors <- colors[1:length(decomp_sig_levels)]
denovo_sig_colors = denovo_sig_colors[denovo_sig_levels %in% unique(tsv$SignatureDenovo)]
decomp_sig_colors = decomp_sig_colors[decomp_sig_levels %in% unique(tsv$SignatureDecomp)]

#reorder points so that Closer/close overlap all others
tsv=tsv; tsv$z="a"
tsv$z[tsv$SVSNV %in% c("CLOSE","CLOSER")] <- "b"
tsv$z[tsv$SNV_ClusterType %in% c("Kataegis","Omikli","DBS","MBS")] <- "b"
tsv <- tsv[order(tsv$z), ]

### SV-SNV clusters SBS6 plot ###
ggplot(data = tsv, aes(x = POS_ABS, y = as.numeric(IMD))) +
  geom_vline(aes(xintercept = POS), data = chr_sizes[1:23,], color="white",size=1.5)+
  scale_x_continuous(breaks=chr_sizes$POS_LABELS,labels=chr_sizes$CHROM,limits = c(-10000000,3105677412),expand = c(0, 0)) +
  ylab("SNV inter-mutation distance (log10)")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position="right")+
  geom_vline(aes(xintercept = POS_ABS,color=TYPE), data = sv,size=1.5,alpha=0.35)+
  scale_color_manual(values=sv_colors)+
  guides(colour = guide_legend(order = 1,nrow=length(unique(sv$TYPE))))+
  new_scale_colour() +
  geom_point(aes(color=MUT6,shape=SVSNV,alpha=z),stat = "identity",size=3.5)+
  scale_alpha_manual(values=c(0.6,1)) +
  scale_shape_manual(values=svsnv_shapes)+
  scale_color_manual(values=snv_colors)+
  scale_y_continuous(trans = "log10",breaks=10^(0:9),labels=formatC(10^(0:9), format="d", big.mark=","))+
  guides(alpha = "none", colour = guide_legend(order = 3,nrow=6), shape = guide_legend(order = 2,nrow=length(svsnv_shapes)))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType)),subtitle = paste0("Closer  ",ocloser," ( ",round(ocloser*1000000/sizecloser,2)," / Mb )",
                                      "  |  Close  ",oclose," ( ",round(oclose*1000000/sizeclose,2)," / Mb )",
                                      "  |  Unclustered  ",ounclustered," ( ",round(ounclustered*1000000/sizeunclustered,2)," / Mb )",
                                      "  |  Breakends  ",nrow(sv))) -> r1

### SNV-SNV clusters plot ### 
ggplot(data = tsv, aes(x = POS_ABS, y = as.numeric(IMD))) +
  geom_vline(aes(xintercept = POS), data = chr_sizes[1:23,], color="white",size=1.5)+
  scale_x_continuous(breaks=chr_sizes$POS_LABELS,labels=chr_sizes$CHROM,limits = c(-10000000,3105677412),expand = c(0, 0)) +
  ylab("SNV inter-mutation distance (log10)")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position="right")+
  geom_vline(aes(xintercept = POS_ABS,color=TYPE), data = sv,size=1.5,alpha=0.35)+
  scale_color_manual(values=sv_colors)+
  guides(colour = guide_legend(order = 1,nrow=length(unique(sv$TYPE))))+
  new_scale_colour() +
  geom_point(aes(color=SNV_ClusterType,shape=SVSNV,alpha=z),stat = "identity",size=3.5)+
  scale_alpha_manual(values=c(0.6,1)) +
  scale_shape_manual(values=svsnv_shapes)+
  scale_color_manual(values=cluster_colors)+
  scale_y_continuous(trans = "log10",breaks=10^(0:9),labels=formatC(10^(0:9), format="d", big.mark=","))+
  guides(alpha = "none", colour = guide_legend(order = 3,nrow=6), shape = guide_legend(order = 2,nrow=length(svsnv_shapes)))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType)),subtitle = paste0("DBS  ",nrow(tsv[which(tsv$SNV_ClusterType=="DBS"),]),
                                                                             "  |  MBS  ",nrow(tsv[which(tsv$SNV_ClusterType=="MBS"),]),
                                                                             "  |  Omikli  ",nrow(tsv[which(tsv$SNV_ClusterType=="Omikli"),])," ( ",sum(tsv$SNV_ClusterType == "Omikli" & tsv$SVSNV == "CLOSER")," Closer / ",
                                                                             sum(tsv$SNV_ClusterType == "Omikli" & tsv$SVSNV == "CLOSE")," Close ) ",
                                                                             "  |  Kataegis  ",nrow(tsv[which(tsv$SNV_ClusterType=="Kataegis"),])," ( ",sum(tsv$SNV_ClusterType == "Kataegis" & tsv$SVSNV == "CLOSER"), " Closer / ",
                                                                             sum(tsv$SNV_ClusterType == "Kataegis" & tsv$SVSNV == "CLOSE")," Close ) ",
                                                                             "  |  Non-clustered  ",nrow(tsv[which(tsv$SNV_ClusterType=="Non-clustered"),])," ( ",sum(tsv$SNV_ClusterType == "Non-clustered" & tsv$SVSNV == "CLOSER")," Closer / ",
                                                                             sum(tsv$SNV_ClusterType == "Non-clustered" & tsv$SVSNV == "CLOSE")," Close ) ")) -> r2

tsv %>% filter(SVSNV=="CLOSER") %>% count(SignatureDenovo) %>% mutate(SignatureDenovo=gsub("-Closer","",SignatureDenovo)) %>% mutate(string = paste0(SignatureDenovo, " - ", n, "  |  ")) %>% summarize(string = paste(string, collapse = "")) %>% mutate(string = str_sub(string, end = -6)) %>% pull(string) -> sub1
tsv %>% filter(SVSNV=="CLOSE") %>% count(SignatureDenovo) %>% mutate(SignatureDenovo=gsub("-Close","",SignatureDenovo)) %>% mutate(string = paste0(SignatureDenovo, " - ", n, "  |  ")) %>% summarize(string = paste(string, collapse = "")) %>% mutate(string = str_sub(string, end = -6)) %>% pull(string) -> sub2
tsv %>% filter(SVSNV=="UNCLUSTERED") %>% count(SignatureDenovo) %>% mutate(SignatureDenovo=gsub("-Unclustered","",SignatureDenovo)) %>% mutate(string = paste0(SignatureDenovo, " - ", n, "  |  ")) %>% summarize(string = paste(string, collapse = "")) %>% mutate(string = str_sub(string, end = -6)) %>% pull(string) -> sub3
### De novo signature plot ### 
ggplot(data = tsv, aes(x = POS_ABS, y = as.numeric(IMD))) +
  geom_vline(aes(xintercept = POS), data = chr_sizes[1:23,], color="white",size=1.5)+
  scale_x_continuous(breaks=chr_sizes$POS_LABELS,labels=chr_sizes$CHROM,limits = c(-10000000,3105677412),expand = c(0, 0)) +
  ylab("SNV inter-mutation distance (log10)")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position="right")+
  geom_vline(aes(xintercept = POS_ABS,color=TYPE), data = sv,size=1.5,alpha=0.35)+
  scale_color_manual(values=sv_colors)+
  guides(colour = guide_legend(order = 1,nrow=length(unique(sv$TYPE))))+
  new_scale_colour() +
  geom_point(aes(color=SignatureDenovo,shape=SVSNV,alpha=z),stat = "identity",size=3.5)+
  scale_alpha_manual(values=c(0.6,1)) +
  scale_shape_manual(values=svsnv_shapes)+
  scale_color_manual(values=denovo_sig_colors)+
  scale_y_continuous(trans = "log10",breaks=10^(0:9),labels=formatC(10^(0:9), format="d", big.mark=","))+
  guides(alpha = "none", colour = guide_legend(order = 3,ncol=1,override.aes=list(shape = legend_shapes)), shape = guide_legend(order = 2,nrow=length(svsnv_shapes)))+
labs(title = paste0(name,"  -  ",unique(tsv$cancerType)),subtitle = paste0("Closer  ",ocloser,"  |  ",sub1,"\n","Close  ",oclose,"  |  ",sub2,"\n","Unclustered  ",ounclustered,"  |  ",sub3)) -> r3

tsv %>% count(SignatureDecomp) %>% mutate(string = paste0(SignatureDecomp, " - ", n, "  |  ")) %>% summarize(string = paste(string, collapse = "")) %>% mutate(string = str_sub(string, end = -6)) %>% pull(string) -> sub4
if (str_count(sub4, "\\|")>11){
  split_string <- str_split(sub4, "  \\|  ")[[1]] 
  count=ceiling((str_count(sub4, "\\|")+1)/2)
  sub4 <- paste0(paste0(split_string[1:count], collapse = "  |  "),"\n",paste0(split_string[(count+1):length(split_string)], collapse = "  |  "))}
### Decomposed signature plot ### 
ggplot(data = tsv, aes(x = POS_ABS, y = as.numeric(IMD))) +
  geom_vline(aes(xintercept = POS), data = chr_sizes[1:23,], color="white",size=1.5)+
  scale_x_continuous(breaks=chr_sizes$POS_LABELS,labels=chr_sizes$CHROM,limits = c(-10000000,3105677412),expand = c(0, 0)) +
  ylab("SNV inter-mutation distance (log10)")+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing = unit(0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position="right")+
  geom_vline(aes(xintercept = POS_ABS,color=TYPE), data = sv,size=1.5,alpha=0.35)+
  scale_color_manual(values=sv_colors)+
  guides(colour = guide_legend(order = 1,nrow=length(unique(sv$TYPE))))+
  new_scale_colour() +
  geom_point(aes(color=SignatureDecomp,shape=SVSNV,alpha=z),stat = "identity",size=3.5)+
  scale_alpha_manual(values=c(0.6,1)) +
  scale_shape_manual(values=svsnv_shapes)+
  scale_color_manual(values=decomp_sig_colors)+
  scale_y_continuous(trans = "log10",breaks=10^(0:9),labels=formatC(10^(0:9), format="d", big.mark=","))+
  guides(alpha = "none", colour = guide_legend(order = 3,ncol=1), shape = guide_legend(order = 2,nrow=length(svsnv_shapes)))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType)),subtitle = sub4) -> r4

tsv %>% select(MUT6,SVSNV) %>% group_by(SVSNV) %>% count(MUT6) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> mut6_plot
ggplot(mut6_plot, aes(fill=MUT6, y=fraction, x=SVSNV)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=snv_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(MUT6), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p1

tsv %>% select(SNV_ClusterType,SVSNV) %>% group_by(SVSNV) %>% count(SNV_ClusterType) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> clusters_plot
ggplot(clusters_plot, aes(fill=SNV_ClusterType, y=fraction, x=SVSNV)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=cluster_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(SNV_ClusterType), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p2

tsv %>% select(SignatureDenovo,SVSNV) %>% group_by(SVSNV) %>% count(SignatureDenovo) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> denovo_sig_plot
ggplot(denovo_sig_plot, aes(fill=SignatureDenovo, y=fraction, x=SVSNV)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=denovo_sig_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, str_extract(denovo_sig_plot$SignatureDenovo, "^[^-]+"), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p3

tsv %>% select(SignatureDecomp,SVSNV) %>% group_by(SVSNV) %>% count(SignatureDecomp) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> decomp_sig_plot
ggplot(decomp_sig_plot, aes(fill=SignatureDecomp, y=fraction, x=SVSNV)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=decomp_sig_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(decomp_sig_plot$SignatureDecomp), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p4

tsv %>% select(MUT6,SNV_ClusterType) %>% group_by(SNV_ClusterType) %>% count(MUT6) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> cluster_mut6_plot
ggplot(cluster_mut6_plot, aes(fill=MUT6, y=fraction, x=SNV_ClusterType)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=snv_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(MUT6), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p5

tsv %>% select(SNV_ClusterType,SVSNV) %>% group_by(SNV_ClusterType) %>% count(SVSNV) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> cluster_svsnv_plot
ggplot(cluster_svsnv_plot, aes(fill=SVSNV, y=fraction, x=SNV_ClusterType)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=svsnv_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(SVSNV), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p6

tsv %>% select(SignatureDenovo,SNV_ClusterType) %>% group_by(SNV_ClusterType) %>% count(SignatureDenovo) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> cluster_denovo_sig_plot
ggplot(cluster_denovo_sig_plot, aes(fill=SignatureDenovo, y=fraction, x=SNV_ClusterType)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=denovo_sig_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(cluster_denovo_sig_plot$SignatureDenovo), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p7

tsv %>% select(SignatureDecomp,SNV_ClusterType) %>% group_by(SNV_ClusterType) %>% count(SignatureDecomp) %>% mutate(total=n %>% sum()) %>% mutate(fraction=100*n/total) %>% select(-total) -> cluster_decomp_sig_plot
ggplot(cluster_decomp_sig_plot, aes(fill=SignatureDecomp, y=fraction, x=SNV_ClusterType)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=decomp_sig_colors)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,25,50,75,100))+
  geom_text(aes(label = ifelse(fraction >5, as.character(cluster_decomp_sig_plot$SignatureDecomp), "")), position = position_stack(vjust = 0.5),size=6) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size = 25),
        plot.title = element_text(hjust = 0.5,size=20,vjust=1),
        panel.background=element_rect(fill="white"),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.spacing = unit(2, "lines"),
        legend.position = "right",
        legend.title = element_blank())+
  ylab("Contribution (%)")+
  guides(fill = guide_legend(ncol=1))+
  labs(title = paste0(name,"  -  ",unique(tsv$cancerType))) -> p8

pdf(file = paste0(name,"_plots.pdf"),width = 23, height = 9.5)
r1; r2; r3; r4
p1 + plot_spacer() + p5 + plot_layout(widths = c(8, 1 ,8))
p2 + plot_spacer() + p6 + plot_layout(widths = c(8, 1 ,8))
p3 + plot_spacer() + p7 + plot_layout(widths = c(8, 1 ,8))
p4 + plot_spacer() + p8 + plot_layout(widths = c(8, 1 ,8))
dev.off()
