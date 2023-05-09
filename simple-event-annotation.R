library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)

#input the name of the vcf via command line; this assumes we input a .sv.vcf.gz file 
#input the name of the sample (cohort, project, sample name)
args=commandArgs(TRUE)
vcf_name=args[1]
sample_name=args[2]
output_vcf_name=paste(sample_name,".sv.ann.vcf",sep="")

# Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))}

#Output annotated vcf
vcf <- readVcf(vcf_name, "hg19")
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
        row.names=c("SIMPLE_TYPE"),
        Number=c("1"),
        Type=c("String"),
        Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))
gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
writeVcf(vcf, output_vcf_name) 

