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

close_bp=params.close_value
closer_bp=params.closer_value
hg19 = file(params.hg19)
fasta_ref=file(params.fasta_ref)
CRG75=file(params.CRG75)
pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()


process make_sv_beds {
       publishDir params.output_folder+"/VCFs/Closer/", mode: 'copy', pattern: '*closer.snv.vcf.gz'
       publishDir params.output_folder+"/VCFs/Close/", mode: 'copy', pattern: '*close.snv.vcf.gz'
       publishDir params.output_folder+"/VCFs/Unclustered/", mode: 'copy', pattern: '*unclustered.snv.vcf.gz'
    
       input:
       set val(sample), file(sv), file(snv) from pairs_list
       file hg19
       file CRG75
       file fasta_ref

       output:
       set file("*closer.snv.vcf.gz") into closer
       set file("*close.snv.vcf.gz") into close
       set file("*unclustered.snv.vcf.gz") into unclustered
       
       shell:
       '''    
       tabix -p vcf !{sv}
       bcftools view -f 'PASS' --regions-file !{CRG75} !{sv} | bcftools sort -Oz > !{sample}.sv.filt.vcf.gz
       bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.filt.vcf.gz > !{sample}.sv.bed
       
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{closer_bp} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.closer.bed
       bedtools slop -i !{sample}.sv.bed -g !{hg19} -b !{close_bp} > !{sample}.cluster.bed
       
       bedtools complement -i !{sample}.cluster.bed -g !{hg19} | sort -k1,1 -k2,2n | bedtools merge > !{sample}.unclustered.bed
       bedtools subtract -a !{sample}.cluster.bed -b !{sample}.closer.bed | sort -k1,1 -k2,2n | bedtools merge > !{sample}.close.bed             
       
       [ -s !{sample}.closer.bed  ] && echo "Closer file not empty" || echo -e '1\t0\t1' >> !{sample}.closer.bed 
       [ -s !{sample}.close.bed  ] && echo "Close file not empty" || echo -e '1\t0\t1' >> !{sample}.close.bed 
       
       tabix -p vcf !{snv}
       bcftools view -f 'PASS' --regions-file !{CRG75} !{snv} | bcftools sort -Oz > !{sample}.filt.vcf.gz
       tabix -p vcf !{sample}.filt.vcf.gz
       
       bcftools view -f PASS --types snps --regions-file !{sample}.closer.bed  !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.!{closer_bp}.closer.snv.vcf.gz
       bcftools view -f PASS --types snps --regions-file !{sample}.close.bed  !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.!{close_bp}.close.snv.vcf.gz
       bcftools view -f PASS --types snps --regions-file !{sample}.unclustered.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Ov > !{sample}.unclustered.snv.vcf.gz
       '''
  }
