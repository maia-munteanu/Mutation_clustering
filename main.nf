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
CRG75 = ile(params.CRG75)
pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()


process get_vcfs {
       publishDir params.output_folder+"/VCFs/Closer_"+closer_bp+"/", mode: 'copy', pattern: '*closer.snv.vcf'
       publishDir params.output_folder+"/VCFs/Close_"+close_bp+"/", mode: 'copy', pattern: '*close.snv.vcf'
       publishDir params.output_folder+"/VCFs/Unclustered/", mode: 'copy', pattern: '*unclustered.snv.vcf'
       publishDir params.output_folder+"/VCFs/Whole/", mode: 'copy', pattern: '*filt.vcf.gz'
    
       input:
       set val sample, path sv, path snv from pairs_list
       path hg19
       path CRG75
       path fasta_ref

       output:
       path "*closer.snv.vcf" into closer
       path "*close.snv.vcf" into close
       path "*unclustered.snv.vcf" into unclustered
       path "*filt.vcf.gz" into others
      
       shell:
       '''    
       Rscript !{baseDir}/simple-event-annotation.R !{sv} !{sample}
       bcftools sort -Oz !{sample}.sv.ann.vcf > !{sample}.sv.ann.vcf.gz
       tabix -p vcf !{sample}.sv.ann.vcf.gz
       bcftools view -f 'PASS' --regions-file !{CRG75} !{sample}.sv.ann.vcf.gz | bcftools sort -Oz > !{sample}.sv.ann.filt.vcf.gz
       bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.ann.filt.vcf.gz > !{sample}.sv.bed
       
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{closer_bp} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.closer.bed
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{close_bp} > !{sample}.cluster.bed
       
       bedtools complement -i !{sample}.cluster.bed -g !{hg19} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.unclustered.bed
       bedtools subtract -a !{sample}.cluster.bed -b !{sample}.closer.bed | sort -k1,1 -k2,2n | bedtools merge > !{sample}.close.bed             
       
       [ -s !{sample}.closer.bed  ] && echo "Closer file not empty" || echo -e '1\t0\t1' >> !{sample}.closer.bed 
       [ -s !{sample}.close.bed  ] && echo "Close file not empty" || echo -e '1\t0\t1' >> !{sample}.close.bed 
       
       tabix -p vcf !{snv}
       bcftools view -f 'PASS' --types snps --regions-file !{CRG75} !{snv} | bcftools sort -Oz > !{sample}.snv.filt.vcf.gz
       tabix -p vcf !{sample}.snv.filt.vcf.gz
       
       bcftools view --regions-file !{sample}.closer.bed  !{sample}.snv.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.closer.snv.vcf
       bcftools view --regions-file !{sample}.close.bed  !{sample}.snv.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.close.snv.vcf
       bcftools view --regions-file !{sample}.unclustered.bed !{sample}.snv.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.unclustered.snv.vcf
       '''
  }
  
process extract96 {
    
    input:
    path "*.vcf" from closer.collect()

    shell:
    '''
    mkdir VCFs && mv *vcf VCFs
    echo "Triggered once after all files complete!"
    python3 !{baseDir}/Extractor.py
   '''     
}
