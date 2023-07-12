#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona

// example run: nextflow run ../nextflow_pipeline/Mutation_clustering/main2.nf --serial_genome /g/strcombio/fsupek_cancer1/SV_clusters_project/Test/hg19.fa.p --chr_sizes /g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome

params.closer_value = 2000
params.close_value = 10000
params.random_window = 1000000
params.random_iter = 5
params.random_batch = 2500
params.svsnv_threshold = 0.2
params.cores = 8
params.sig_cores = 32
params.minsig = 1
params.maxsig = 1
params.assembly = "hg19"
params.sigproassembly = "GRCh37"


params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input2.tsv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/Main2Results"
params.mappability = "/home/mmunteanu/reference/CRG75_nochr.bed"
params.reference = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.serial_genome = "/g/strcombio/fsupek_cancer1/SV_clusters_project/nextflow_analysis/work/fe/de893a2682e0e596c0e93511326ac1/hg19.fa.p"
params.chr_sizes = "/g/strcombio/fsupek_cancer1/SV_clusters_project/nextflow_analysis/work/74/ef2caf9be07977563e32c228dc4bab/hg19.genome"
params.vcfanno_conf = "/g/strcombio/fsupek_cancer1/SV_clusters_project/vcfanno/vcfanno.conf"

input_file = file(params.input_file)
reference = file(params.reference)
mappability = file(params.mappability)
vcfanno_conf = file(params.vcfanno_conf)

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

sv_list = Channel.fromPath(input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.sv), file(row.linx) ] }

process parse_svs {
       tag { sample }
       input:
       tuple val(sample), file(sv), file(linx) from sv_list
       path mappability
       path chr from chr_sizes
       
       output:
       tuple val(sample), file("${sample}.sv_snv.ann.bed"), optional: true into(svs_exist, filter_by_sv_snv) 
       tuple val(sample), file("${sample}.sv.ann.tsv"), optional: true into annotate_with_sv_info
       tuple val(sample), env(filter) into sv_filter
      
       shell:
       '''  
       if [[ $(zgrep -v "^#" !{sv} | wc -l) -gt 0 ]]
       then
              svname=$(bcftools query -l !{sv} | sed -n 2p)
              Rscript !{baseDir}/simple-event-annotation.R !{sv} !{sample}
              bcftools sort -Oz !{sample}.sv.ann.vcf > !{sample}.sv.ann.vcf.gz
              tabix -p vcf !{sample}.sv.ann.vcf.gz
              bcftools view -s $svname -f 'PASS' --regions-file !{mappability} !{sample}.sv.ann.vcf.gz | bcftools sort -Oz > !{sample}.sv.ann.filt.vcf.gz
              
              if [[ $(zgrep -v "^#" !{sample}.sv.ann.filt.vcf.gz | wc -l) -gt 0 ]] 
              then             
                     bcftools query -f '%CHROM\t%POS\t%ID\t%QUAL\t%SVLEN\t%SIMPLE_TYPE[\t%PURPLE_AF][\t%PURPLE_CN]\n' !{sample}.sv.ann.filt.vcf.gz > !{sample}.sv.ann.tsv
                     Rscript !{baseDir}/linx_annotation.R !{sample}.sv.ann.tsv !{linx}
                     
                     bcftools query -f '%CHROM\t%POS\t%POS\n' !{sample}.sv.ann.filt.vcf.gz > sv.bed
                     bedtools slop -i sv.bed -g !{chr} -b !{params.closer_value} | sort -k1,1 -k2,2n | bedtools merge > closer.bed
                     bedtools slop -i sv.bed -g !{chr} -b !{params.close_value} > cluster.bed
                     bedtools complement -i cluster.bed -g !{chr} | sort -k1,1 -k2,2n | bedtools merge > unclustered.bed
                     bedtools subtract -a cluster.bed -b closer.bed | sort -k1,1 -k2,2n | bedtools merge > close.bed     
                     awk -v OFS='\t' '{print $1,$2,$3,"CLOSER"}' closer.bed > closer.ann.bed
                     awk -v OFS='\t' '{print $1,$2,$3,"CLOSE"}' close.bed > close.ann.bed
                     awk -v OFS='\t' '{print $1,$2,$3,"UNCLUSTERED"}' unclustered.bed > unclustered.ann.bed
                     cat *ann.bed | sort -k 1,1 -k2,2n > !{sample}.sv_snv.ann.bed 
                     filter="SV-FILTER-PASS"
              else
                     filter="SV-FILTER-2"
              fi
       else 
              filter="SV-FILTER-1"  
       fi        
       '''
}

snv_list = Channel.fromPath(input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.join(svs_exist).map{ sample, snv, bed -> [ sample, snv ] }

process parse_snvs {
       tag { sample }
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
       tag { sample }
       errorStrategy 'retry'
       memory { 30.GB * task.attempt }
       input:
       path serial from serial_genome
       tuple val(sample), file(snv2rand) from snvs_to_randomise
      
       output:
       tuple val(sample), file("${sample}.snv.filt.vcf.gz"), file("${sample}.snv.filt.random*.vcf") into randomised_vcf
       tuple val(sample), file("${sample}.snv.filt.random*.tsv") into randomised_tsv
       
       shell:
       '''
       bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT{0}\t1\t!{sample}\n' !{snv2rand} > !{sample}.snv.filt.tsv
       randommut -M randomize -g !{serial} -m !{sample}.snv.filt.tsv -o !{sample}.snv.filt.random.R!{params.random_iter}.tsv -t !{params.random_iter} -w !{params.random_window} -b !{params.random_batch}
       Rscript !{baseDir}/tsv_to_vcf.R !{sample} !{params.random_iter} !{sample}.snv.filt.tsv !{sample}.snv.filt.random.R!{params.random_iter}.tsv !{sample}.snv.filt.random.R!{params.random_iter}.vcf
       '''
}

sv_snv = randomised_vcf.join(filter_by_sv_snv)

process get_sv_snv_clusters {
       tag { sample }
       errorStrategy 'retry'
       memory { 10.GB * task.attempt }
       input:
       tuple val(sample), file(ovcf), file(rvcf), file(bed) from sv_snv
      
       output:
       tuple val(sample), file("${sample}.snv.filt.svsnv.vcf.gz"), optional: true into annotate_snvs
       tuple file("${sample}.snv.closer.vcf"), file("${sample}.snv.close.vcf"), file("${sample}.snv.unclustered.vcf"), optional: true into to_count
       tuple val(sample), env(filter), env(ratio), env(rcloser), env(rclose), env(runclustered), env(ocloser), env(oclose), env(ounclustered), env(sizecloser), env(sizeclose), env(sizeunclustered) into(snv_filter, sample_info)

       shell:
       '''
       sizecloser=$(grep -w CLOSER !{sample}.sv_snv.ann.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{printf "%.0f", SUM}')
       sizeclose=$(grep -w CLOSE !{sample}.sv_snv.ann.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{printf "%.0f", SUM}')
       sizeunclustered=$(grep -w UNCLUSTERED !{sample}.sv_snv.ann.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=\$3-\$2 }END{printf "%.0f", SUM}')
       
       bgzip !{sample}.sv_snv.ann.bed
       tabix -p bed !{sample}.sv_snv.ann.bed.gz
       echo '[[annotation]] \n file=\"!{sample}.sv_snv.ann.bed.gz\" \n names=[\"SVSNV\"] \n columns=[4] \n ops=[\"self\"]' >> !{sample}.conf
       vcfanno_linux64 !{sample}.conf !{rvcf} > !{sample}.snv.filt.random.R!{params.random_iter}.svsnv.vcf
       vcfanno_linux64 !{sample}.conf !{ovcf} > !{sample}.snv.filt.svsnv.vcf
       bgzip !{sample}.snv.filt.random.R!{params.random_iter}.svsnv.vcf; bgzip !{sample}.snv.filt.svsnv.vcf
       
       ocloser=$(zgrep -w SVSNV=CLOSER !{sample}.snv.filt.svsnv.vcf.gz | wc -l); oclose=$(zgrep -w SVSNV=CLOSE !{sample}.snv.filt.svsnv.vcf.gz | wc -l); ounclustered=$(zgrep -w SVSNV=UNCLUSTERED !{sample}.snv.filt.svsnv.vcf.gz | wc -l)
       rcloser=$(zgrep -w SVSNV=CLOSER !{sample}.snv.filt.random.R!{params.random_iter}.svsnv.vcf.gz | wc -l); rclose=$(zgrep -w SVSNV=CLOSE !{sample}.snv.filt.random.R!{params.random_iter}.svsnv.vcf.gz | wc -l); runclustered=$(zgrep -w SVSNV=UNCLUSTERED !{sample}.snv.filt.random.R!{params.random_iter}.svsnv.vcf.gz | wc -l)
       echo "Observed closer n=$ocloser"; echo "Observed close n=$oclose"; echo "Randomised closer n=$rcloser"; echo "Randomised close n=$rclose"

       if [[ $ocloser -gt 0 && $oclose -gt 0 ]]
       then
             ratio=$(echo "scale=3; ($rcloser+$rclose)/(($ocloser+$oclose)*!{params.random_iter})" | bc); echo "Ratio is $ratio"
             if [[ $(echo "$ratio<=!{params.svsnv_threshold}" | bc) -eq 1 ]]
             then 
                   echo "Sample passes filters 1. SV-SNV are present and 2. below !{params.svsnv_threshold} randomised clusters."
                   bcftools filter -i 'INFO/SVSNV=="CLOSER"' !{sample}.snv.filt.svsnv.vcf.gz | bcftools sort -Ov > !{sample}.snv.closer.vcf
                   bcftools filter -i 'INFO/SVSNV=="CLOSE"' !{sample}.snv.filt.svsnv.vcf.gz | bcftools sort -Ov > !{sample}.snv.close.vcf
                   bcftools filter -i 'INFO/SVSNV=="UNCLUSTERED"' !{sample}.snv.filt.svsnv.vcf.gz | bcftools sort -Ov > !{sample}.snv.unclustered.vcf
                   filter="SNV-FILTER-PASS"
             else
                   echo "Sample passes filter 1. SV-SNV are present but fails filter 2. below !{params.svsnv_threshold} randomised clusters."; rm !{sample}.snv.filt.svsnv.vcf.gz
                   filter="SNV-FILTER-2"
                   
             fi
       else
             echo "Sample does not have SV-SNV clusters"; rm !{sample}.snv.filt.svsnv.vcf.gz; 
             filter="SNV-FILTER-1"; ratio=0
       fi  
       echo $sizecloser 
       echo $sizeclose 
       echo $sizeunclustered
       '''
}

process get_snv_clusters {
       tag { sample }
       publishDir params.output_folder+"/Plots/SNV-SNV", mode: 'move', pattern: '*_plot.pdf'
    
       input:
       tuple val(sample), file(tsv) from randomised_tsv 
      
       output:
       tuple val(sample), file("${sample}.snv.clusters.tsv") into snv_clusters 
       file("${sample}_plot.pdf")
       
       shell:
       '''
       clustmut distance -i . --glob !{tsv} -o !{sample} -Vv
       Rscript !{baseDir}/vranges_to_tsv.R !{sample} !{sample}_distance_VRanges.rds
       '''
}

process count_mutations {
    publishDir params.output_folder+"/Counts/", mode: 'copy', pattern: '*.all'
    
    input:
    path "*" from to_count.collect()

    output:
    path("*.all") into counts_all

    shell:
    '''
    mkdir closer_VCFs && mv *snv.closer.vcf closer_VCFs
    mkdir close_VCFs && mv *snv.close.vcf close_VCFs
    mkdir unclustered_VCFs && mv *snv.unclustered.vcf unclustered_VCFs
    
    python3 !{baseDir}/MatrixGenerator.py "Closer" !{params.sigproassembly} "./closer_VCFs/"
    python3 !{baseDir}/MatrixGenerator.py "Close" !{params.sigproassembly} "./close_VCFs/"
    python3 !{baseDir}/MatrixGenerator.py "Unclustered" !{params.sigproassembly} "./unclustered_VCFs/"
    
    cp ./closer_VCFs/output/SBS/Closer.SBS96.all ./Closer.all
    cp ./close_VCFs/output/SBS/Close.SBS96.all ./Close.all
    cp ./unclustered_VCFs/output/SBS/Unclustered.SBS96.all ./Unclustered.all
   '''     
}

counts = counts_all.flatten().map { file -> [file.baseName, file] }.view()

process get_signatures {
    tag { name }
    cpus = params.sig_cores

    publishDir params.output_folder+"/Signatures/", mode: 'move', pattern: "${name}"

    input:
    tuple val(name), path(count) from counts

    output:
    tuple path("*denovo.txt"), path("*decomp.txt") into probabilities
    path("${name}")

    shell:
    '''
    python3 !{baseDir}/SignatureExtractor.py "./Signatures" !{count} !{params.sigproassembly} !{params.minsig} !{params.maxsig} !{params.sig_cores}
    cp ./Signatures/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/De_Novo_Mutation_Probabilities_refit.txt ./!{name}_denovo.txt
    cp ./Signatures/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt ./!{name}_decomp.txt
    mv ./Signatures ./!{name}
    '''
}

snv_to_annotate = annotate_snvs.join(snv_clusters).join(annotate_with_sv_info).join(sample_info)

closer_decomp=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Closer_mutations_output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt")
closer_denovo=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Closer_mutations_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/De_Novo_Mutation_Probabilities_refit.txt")
close_decomp=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Close_mutations_output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt")
close_denovo=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Close_mutations_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/De_Novo_Mutation_Probabilities_refit.txt")
unclustered_decomp=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Unclustered_mutations_output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/Decomposed_Mutation_Probabilities.txt")
unclustered_denovo=file("/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Unclustered_mutations_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/De_Novo_Mutation_Probabilities_refit.txt")

process snv_annotation {
       tag { sample }
       cpus = params.cores
       errorStrategy 'retry'
       memory { 10.GB * task.attempt }

       publishDir params.output_folder+"/Plots/SNVs", mode: 'move', pattern: '*_plots.pdf'
       publishDir params.output_folder+"/Annotated-SNVs", mode: 'move', pattern: '*_annotated.tsv'
    
       input:
       tuple val(sample), file(vcf), file(snvsnv), file(sv), val(filter), val(ratio), val(rcloser), val(rclose), val(runclustered), val(ocloser), val(oclose), val(ounclustered), val(sizecloser), val(sizeclose), val(sizeunclustered) from snv_to_annotate 
       path "*" from probabilities.collect()
       path input_file  
       path chr from chr_sizes
       path vcfanno_conf

       output:
       tuple val(sample), file("${sample}_annotated.tsv"), file("${sample}_plots.pdf") 
       
       shell:
       '''
       vcfanno_linux64  !{vcfanno_conf} !{vcf} > !{sample}.snv.filt.svsnv.ann.vcf
       bgzip !{sample}.snv.filt.svsnv.ann.vcf
       tabix -p vcf !{sample}.snv.filt.svsnv.ann.vcf.gz
       vcf2tsv -n NA !{sample}.snv.filt.svsnv.ann.vcf.gz > !{sample}.snv.filt.svsnv.ann.tsv
       #Rscript !{baseDir}/snv_annotation.R !{params.cores} !{params.closer_value} !{params.close_value} !{input_file} !{params.assembly} !{chr} !{sample} !{sample}.snv.filt.svsnv.ann.tsv !{snvsnv} !{sv} !{ratio} !{ocloser} !{oclose} !{ounclustered} !{sizecloser} !{sizeclose} !{sizeunclustered} Closer_decomp.txt Closer_denovo.txt Close_decomp.txt Close_denovo.txt Unclustered_decomp.txt Unclustered_denovo.txt
       Rscript !{baseDir}/snv_annotation.R !{params.cores} !{params.closer_value} !{params.close_value} !{input_file} !{params.assembly} !{chr} !{sample} !{sample}.snv.filt.svsnv.ann.tsv !{snvsnv} !{sv} !{ratio} !{ocloser} !{oclose} !{ounclustered} !{sizecloser} !{sizeclose} !{sizeunclustered} !{closer_decomp} !{closer_denovo} !{close_decomp} !{close_denovo} !{unclustered_decomp} !{unclustered_denovo}

       '''
}

sv_filter.flatten()
               .collate( 2 )
               .map { it.join("\t") + "\n" } 
               .collectFile(name: 'SV-filter.tsv', storeDir: params.output_folder+"/Filtering_outcome/")

snv_filter.flatten()
               .collate( 12 )
               .map { it.join("\t") + "\n" } 
               .collectFile(name: 'SNV-filter.tsv', storeDir: params.output_folder+"/Filtering_outcome/")


    


