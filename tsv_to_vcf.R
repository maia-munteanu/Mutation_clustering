library(data.table)
args=commandArgs(TRUE)

sample=args[1]
niter=args[2]
otsv=args[3]
rtsv=args[4]
output=args[5]

chromosomes=c(1:22,"X","Y")

original_tsv=read.table(file = otsv, sep = '\t', header = FALSE)
colnames(original_tsv)=c("chr","start","end","ref","alt","strand","sample")
random_tsv=read.table(file = rtsv, sep = '\t', header = TRUE)

vcf=data.frame(CHROM = c(), POS = c(), ID=c(), REF = c(), ALT = c(),  QUAL=c(), FILTER=c(), INFO=c())
if(identical(original_tsv[,c(1,3,4,5)],random_tsv[,c(1,3,6,7)])){ #check that chr, pos, ref and alt fields for the original variants are the same
      for (i in 1:niter){vcf=rbind(vcf, data.frame(CHROM = random_tsv$chr, POS = random_tsv[,paste0("R",i)], ID="rsX", REF = random_tsv$ref, ALT = random_tsv$alt,  QUAL=".", FILTER="PASS",INFO=paste0("R",i)) )}
      vcf$CHROM=factor(vcf$CHROM,levels=chromosomes[chromosomes %in% unique(vcf$CHROM)])
      vcf=vcf[order(vcf$CHROM,vcf$POS ),]
      colnames(vcf)[which(grepl("CHROM", colnames(vcf)))] = "#CHROM"
      write.table(c("##fileformat=VCFv4.2"), output, row.names = F,col.names = F,quote = F)
      fwrite(vcf, output, row.names = F,col.names = T,quote = F,append = T,sep="\t")
}else{write("The observed variant information found in the randomised tsv file do not match the original variant information. Please check the output of randommut.", stderr())}
