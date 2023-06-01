#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona

// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main2.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p --chr_sizes /g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome

params.closer_value = 2000
params.close_value = 10000
params.random_window = 1000000
params.random_iter = 5
params.random_batch = 2500
params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input2.csv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Main2Results"
params.mappability = "/home/mmunteanu/reference/CRG75_nochr.bed"
params.reference = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.assembly = "hg19"
params.serial_genome = null
params.chr_sizes = null

reference = file(params.reference)
mappability = file(params.mappability)

if (params.serial_genome){
      serial_genome = file(params.serial_genome)
}else{      
    process serialize_genome {
        input: 
        path "${params.assembly}.fa" from reference

        output:
        path "${params.assembly}.fa.p" into serial_genome

        """
        randommut -M serialize -g ${params.assembly}.fa -a ${params.assembly}
        """
    }
}
    
if (params.chr_sizes){
      chr_sizes = file(params.chr_sizes)
}else{      
    process get_chr_sizes {
        input: 
        path "${params.assembly}.fa" from reference

        output:
        path "${params.assembly}.genome" into chr_sizes

        """
        samtools faidx ${params.assembly}.fa
        head -n 24 ${params.assembly}.fa.fai | cut -f1,2  > ${params.assembly}.genome
        """
    }
}

sv_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv) ] }.view()

process parse_svs {
       input:
       tuple val(sample), file(sv) from sv_list
       path mappability
       path chr_sizes
       
       output:
       tuple val(sample), file("${sample}.sv_snv.ann.bed"), optional: true into(svs_exist, filter_by_sv_snv) 
       tuple val(sample), file("${sample}.sv.ann.txt"), optional: true into annotate_with_sv_info
      
       shell:
       '''  
       if [ $(zgrep -v "^#" !{sv} | wc -l) -gt 0 ]
       then
              svname=$(bcftools query -l !{sv} | sed -n 2p)
              Rscript !{baseDir}/simple-event-annotation.R !{sv} !{sample}
              bcftools sort -Oz !{sample}.sv.ann.vcf > !{sample}.sv.ann.vcf.gz
              tabix -p vcf !{sample}.sv.ann.vcf.gz
              bcftools view -s $svname -f 'PASS' --regions-file !{mappability} !{sample}.sv.ann.vcf.gz | bcftools sort -Oz > !{sample}.sv.ann.filt.vcf.gz
       fi        
           
       if [ -f !{sample}.sv.ann.filt.vcf.gz ] && [ $(zgrep -v "^#" !{sample}.sv.ann.filt.vcf.gz | wc -l) -gt 0 ] 
       then            
             bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.ann.filt.vcf.gz > sv.bed
             bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SIMPLE_TYPE\n' !{sample}.sv.ann.filt.vcf.gz > !{sample}.sv.ann.txt

             bedtools slop -i sv.bed -g !{chr_sizes} -b !{params.closer_value} | sort -k1,1 -k2,2n | bedtools merge > closer.bed
             bedtools slop -i sv.bed -g !{chr_sizes} -b !{params.close_value} > cluster.bed
             bedtools complement -i cluster.bed -g !{chr_sizes} | sort -k1,1 -k2,2n | bedtools merge > unclustered.bed
             bedtools subtract -a cluster.bed -b closer.bed | sort -k1,1 -k2,2n | bedtools merge > close.bed     
             awk -v OFS='\t' '{print $1,$2,$3,"CLOSER"}' closer.bed > closer.ann.bed
             awk -v OFS='\t' '{print $1,$2,$3,"CLOSE"}' close.bed > close.ann.bed
             awk -v OFS='\t' '{print $1,$2,$3,"UNCLUSTERED"}' unclustered.bed > unclustered.ann.bed
             cat *ann.bed | sort -k 1,1 -k2,2n > !{sample}.sv_snv.ann.bed 
       fi
       '''
}

snv_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.join(svs_exist).view()

process parse_snvs {
       input:
       tuple val(sample), file(snv) from snv_list
       path mappability
       
       output:
       tuple val(sample), file("${sample}.snv.filt.vcf.gz") into snvs_to_randomise
      
       shell:
       '''  
       snvname=$(bcftools query -l !{snv} | sed -n 2p)
       tabix -p vcf !{snv}
       bcftools view -s $snvname -f 'PASS' --types snps --regions-file !{mappability} !{snv} | bcftools sort -Oz > !{sample}.snv.filt.vcf.gz
       tabix -p vcf !{sample}.snv.filt.vcf.gz       
       '''
}

process randomise_snvs {
errorStrategy 'retry'
       memory { 30.GB * task.attempt }
      
       input:
       serial_genome
       tuple val(sample), file(snv2rand) from snvs_to_randomise
      
       output:
       tuple val(sample), file("${sample}.snv.filt.vcf.gz"), file("${sample}.snv.filt.random*.vcf") into randomised_vcf
       tuple val(sample), file("${sample}.snv.filt.random*.tsv") into randomised_tsv

       
       shell:
       '''
       bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT{0}\t1\t!{sample}\n' !{snv2rand} > !{sample}.snv.filt.tsv
       randommut -M randomize -g !{serial_genome} -m !{sample}.snv.filt.tsv -o !{sample}.snv.filt.random.R!{params.random_iter}.tsv -t !{params.random_iter} -w !{params.random_window} -b !{params.random_batch}
       Rscript !{baseDir}/tsv_to_vcf.R !{sample} !{params.random_iter} !{sample}.snv.filt.tsv !{sample}.snv.filt.random.R!{params.random_iter}.tsv !{sample}.snv.filt.random.R!{params.random_iter}.vcf

       '''
}

process get_snv_clusters {
      
       input:
       tuple val(sample), file(tsv) from randomised_tsv 
      
       output:
       tuple val(sample), file("${sample}.snv.clusters.tsv") into snv_clusters 
       
       shell:
       '''
       clustmut distance -i . --glob !{tsv} -o !{sample} -Vv
       Rscript !{baseDir}/vranges_to_tsv.R !{sample} !{sample}_distance_VRanges.rds
       '''
}

//tuple val(sample), file("${sample}.sv_snv.ann.bed") into filter_by_sv_snv 

sv_snv = randomised_vcf.join(filter_by_sv_snv).view()

process get_sv_snv_clusters {
       input:
       tuple val(sample), file(ovcf), file(rvcf), file(bed) from sv_snv      
       
       output:
       tuple val(sample), file("${sample}.sv_snv.ann.bed.gz.tbi")
             
       shell:
       '''
       bgzip !{sample}.sv_snv.ann.bed
       tabix -p bed !{sample}.sv_snv.ann.bed.gz
       '''
      
       
       
       
}



