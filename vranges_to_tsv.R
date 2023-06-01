args=commandArgs(TRUE)

sample=args[1]
vr=args[2]

rds=readRDS(file=vr)

df = data.frame(
     sample=sampleNames(rds),
     chr=seqnames(rds),
     pos=start(rds),
     ref=ref(rds),
     alt=alt(rds),
     sbs96=rds$ctx,
     fdr=rds$fdr,
     distance=rds$dist,
     exp_distance=rds$exp_dist,
     cluster_type=rds$event_type,
     cluster_size=rds$event_muts,
     cluster_id=rds$event_rid, stringsAsFactors = F)
readr::write_tsv(df,file = paste0(sample,".snv.clusters.tsv"))

