#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona


params.closer_value = 2000
params.close_value = 10000
params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input.csv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/New_pipe_results"
params.fasta_ref = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.hg19 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome"
params.CRG75 = "/home/mmunteanu/reference/CRG75_nochr.bed"

close_bp = params.close_value
closer_bp = params.closer_value
fasta_ref = file(params.fasta_ref)
hg19 = file(params.hg19)
CRG75 = file(params.CRG75)
pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()

process parse_vcfs {
       publishDir params.output_folder+"/VCFs/Closer_"+closer_bp+"/", mode: 'copy', pattern: '*closer.snv.vcf'
       publishDir params.output_folder+"/VCFs/Close_"+close_bp+"/", mode: 'copy', pattern: '*close.snv.vcf'
       publishDir params.output_folder+"/VCFs/Unclustered/", mode: 'copy', pattern: '*unclustered.snv.vcf'
       publishDir params.output_folder+"/VCFs/Whole/", mode: 'copy', pattern: '*filt.vcf.gz'
    
       input:
       tuple val(sample), path(sv), path(snv) from pairs_list
       path hg19
       path CRG75
       path fasta_ref

       output:
       path "*closer.snv.vcf" into closer
       path "*close.snv.vcf" into close
       path "*unclustered.snv.vcf" into unclustered
       path "*.snv.filt.vcf.gz" into whole
      
       shell:
       '''  
       svname=$(bcftools query -l !{sv} | sed -n 2p)
       snvname=$(bcftools query -l !{snv} | sed -n 2p)
       
       Rscript !{baseDir}/simple-event-annotation.R !{sv} !{sample}
       bcftools sort -Oz !{sample}.sv.ann.vcf > !{sample}.sv.ann.vcf.gz
       tabix -p vcf !{sample}.sv.ann.vcf.gz
       bcftools view -s $svname -f 'PASS' --regions-file !{CRG75} !{sample}.sv.ann.vcf.gz | bcftools sort -Oz > !{sample}.sv.ann.filt.vcf.gz
       bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.ann.filt.vcf.gz > !{sample}.sv.bed
       
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{closer_bp} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.closer.bed
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{close_bp} > !{sample}.cluster.bed
       bedtools complement -i !{sample}.cluster.bed -g !{hg19} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.unclustered.bed
       bedtools subtract -a !{sample}.cluster.bed -b !{sample}.closer.bed | sort -k1,1 -k2,2n | bedtools merge > !{sample}.close.bed             
       [ -s !{sample}.closer.bed  ] && echo "Closer file not empty" || echo -e '1\t0\t1' >> !{sample}.closer.bed 
       [ -s !{sample}.close.bed  ] && echo "Close file not empty" || echo -e '1\t0\t1' >> !{sample}.close.bed 
       
       tabix -p vcf !{snv}
       bcftools view -s $snvname -f 'PASS' --types snps --regions-file !{CRG75} !{snv} | bcftools sort -Oz > !{sample}.snv.filt.vcf.gz
       tabix -p vcf !{sample}.snv.filt.vcf.gz
       
       bcftools view --regions-file !{sample}.closer.bed  !{sample}.snv.filt.vcf.gz | bcftools norm -d none -f !{fasta_ref} | bcftools sort -Ov > !{sample}.closer.snv.vcf
       bcftools view --regions-file !{sample}.close.bed  !{sample}.snv.filt.vcf.gz |  bcftools norm -d none -f !{fasta_ref} | bcftools sort -Ov > !{sample}.close.snv.vcf
       bcftools view --regions-file !{sample}.unclustered.bed !{sample}.snv.filt.vcf.gz |  bcftools norm -d none -f !{fasta_ref} | bcftools sort -Ov > !{sample}.unclustered.snv.vcf
       '''
  }
  
process get_clusters {
    input:
    path "*" from whole.collect()
    
    shell:
    '''  
    gunzip --force *.vcf.gz
    mkdir VCFs && mv *snv.filt.vcf VCFs
    python3 !{baseDir}/ClusterClassifier.py "clusters" "GRCh37" "./VCFs/"
    ''' 
}
 
process count_mutations {
    publishDir params.output_folder+"/Counts/", mode: 'copy', pattern: '*.all'
    
    input:
    path "*" from closer.collect()
    path "*" from close.collect()
    path "*" from unclustered.collect()
    
    output:
    path "*.all" into counts

    shell:
    '''
    mkdir closer_VCFs && mv *closer.snv.vcf closer_VCFs
    mkdir close_VCFs && mv *close.snv.vcf close_VCFs
    mkdir unclustered_VCFs && mv *unclustered.snv.vcf unclustered_VCFs
    
    python3 !{baseDir}/MatrixGenerator.py "closer" "GRCh37" "./closer_VCFs/"
    python3 !{baseDir}/MatrixGenerator.py "close" "GRCh37" "./close_VCFs/"
    python3 !{baseDir}/MatrixGenerator.py "unclustered" "GRCh37" "./unclustered_VCFs/"
    
    cp ./closer_VCFs/output/SBS/closer.SBS96.all ./
    cp ./close_VCFs/output/SBS/close.SBS96.all ./
    cp ./unclustered_VCFs/output/SBS/unclustered.SBS96.all ./
   '''     
}

process get_signatures {
    input:
    path "*" from counts
    
    shell:
    '''
    mkdir Unclustered && mv closer.SBS96.all Closer
    mkdir Close && mv close.SBS96.all Close
    mkdir Unclustered && mv unclustered.SBS96.all Unclustered
    
    python3 !{baseDir}/SignatureExtractor.py "./Closer/Signatures" "./Closer/closer.SBS96.all" "GRCh37" 1 5
    python3 !{baseDir}/SignatureExtractor.py "./Close/Signatures" "./Close/close.SBS96.all" "GRCh37" 1 5
    python3 !{baseDir}/SignatureExtractor.py "./Unclustered/Signatures" "./Unclustered/unclustered.SBS96.all" "GRCh37" 1 5
    '''
}
