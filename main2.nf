#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona

// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main2.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p --chr_sizes /g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome

params.closer_value = 2000
params.close_value = 10000
params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input.csv"
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

samples_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.snv) ] }.view()
                   
                   
process parse_vcfs {
       input:
       tuple val(sample), file(sv), file(snv) from samples_list
       path mappability
       path chr_sizes
       
       output:
       set val(sample), file("${sample}.snv.filt.vcf.gz") into snvs_to_randomise
         
       shell:
       '''  
       svname=$(bcftools query -l !{sv} | sed -n 2p)
       snvname=$(bcftools query -l !{snv} | sed -n 2p)
       
       Rscript !{baseDir}/simple-event-annotation.R !{sv} !{sample}
       bcftools sort -Oz !{sample}.sv.ann.vcf > !{sample}.sv.ann.vcf.gz
       tabix -p vcf !{sample}.sv.ann.vcf.gz
       bcftools view -s $svname -f 'PASS' --regions-file !{mappability} !{sample}.sv.ann.vcf.gz | bcftools sort -Oz > !{sample}.sv.ann.filt.vcf.gz
       bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.ann.filt.vcf.gz > sv.bed
       bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SIMPLE_TYPE\n' !{sample}.sv.ann.filt.vcf.gz > !{sample}.sv.ann.txt
       
       bedtools slop -i sv.bed -g !{chr_sizes} -b !{params.closer_value} | sort -k1,1 -k2,2n | bedtools merge > closer.bed
       bedtools slop -i sv.bed -g !{chr_sizes} -b !{params.close_value} > cluster.bed
       bedtools complement -i cluster.bed -g !{chr_sizes} | sort -k1,1 -k2,2n | bedtools merge > unclustered.bed
       bedtools subtract -a cluster.bed -b closer.bed | sort -k1,1 -k2,2n | bedtools merge > close.bed     
       
       [ -s closer.bed  ] && echo "Closer file not empty" || echo -e '1\t0\t1' >> closer.bed 
       [ -s close.bed  ] && echo "Close file not empty" || echo -e '1\t0\t1' >> close.bed 
       
       tabix -p vcf !{snv}
       bcftools view -s $snvname -f 'PASS' --types snps --regions-file !{mappability} !{snv} | bcftools sort -Oz > !{sample}.snv.filt.vcf.gz
       tabix -p vcf !{sample}.snv.filt.vcf.gz       
       '''
       
}


process randomise_snvs {
       errorStrategy 'retry'
       memory { 30GB * task.attempt }
      
       input:
       serial_genome
       set val(sample), file(snv2rand) from snvs_to_randomise
       
       shell:
       '''
       bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT{0}\t1\tsampleA\n' !{snv2rand} > !{sample}.snv.tsv
       randommut -M randomize -g !{serial_genome} -m !{sample}.snv.tsv -o !{sample}.random.snv.tsv -t 1 -w 1000000 -b 2500
       '''

}
