#getting variables from bash
spikein_genome<-Sys.getenv("SPIKEIN_GENOME")
data_info_file<-Sys.getenv("DATA_INFO_FILE")
aligned_bam_dir<-Sys.getenv("ALIGNED_BAM_DIR")
thread_count<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
spikein_gtf_file<-Sys.getenv("SPIKEIN_GTF_FILE")
fcounts_deseq_dir<-Sys.getenv("FC_DSQ_DIR")
refr_group_name<-Sys.getenv("REFR_GROUP_NAME")
se_or_pe<-Sys.getenv("SE_or_PE")
seq_stranded<-Sys.getenv("SEQ_STRANDED")
norm_factor_file<-Sys.getenv("NORM_FACTOR_FILE")
# factor_name<-Sys.getenv("FACTOR_NAME")

#set up sampleSheet for spike_genome counting
library(data.table)
pre_samplesheet<-fread(data_info_file)
pre_samplesheet
samplesheet<-unique(pre_samplesheet[,.(SampleID,Condition,Replicate)])
samplesheet[,bamReads:=paste0(aligned_bam_dir,SampleID,"_",spikein_genome,"_Aligned.sortedByCoord.out.bam"),by=SampleID]

#samplesheet[,Factor:=factor_name]

# setting variables for featureCounts counting
if(se_or_pe=="PE"){
	paired_setting=T
} else if (se_or_pe=="SE"){
	paired_setting=F
}

if(seq_stranded=="Yes"){
	strand_setting<-2
} else if(seq_stranded=="No"){
	strand_setting<-0
}

cat("Starting featureCounts\n")
library(Rsubread)
fcounts<-featureCounts(files=samplesheet$bamReads, 
	annot.ext=spikein_gtf_file, isPairedEnd=paired_setting, nthreads=thread_count,
	isGTFAnnotationFile=T, strandSpecific=strand_setting)

saveRDS(fcounts,paste0(fcounts_deseq_dir,"spikein_fcounts_output.rds"))

# adjust colnames for DESeq2 setup
samplesheet[,fcounts_colname:=basename(bamReads),by=SampleID]
colnames(fcounts$counts)<- sapply(colnames(fcounts$counts),
                        function(x){samplesheet[fcounts_colname==x,SampleID]})

# making colData info for DESeq2
colData<-as.data.frame(samplesheet[,.(Condition)])
rownames(colData)<-samplesheet$SampleID
colData$Condition<-factor(colData$Condition)
colData$Condition<-relevel(colData$Condition,ref=refr_group_name)
colData

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData=fcounts$counts,colData=colData,design=~Condition)

# run DESeq
dds<-DESeq(dds)
res<-results(dds)
summary(res)

# get 1/sizeFactors = normFactors
norm_factors<-1/sizeFactors(dds)
norm_factor_dt<-data.table(SampleID=names(norm_factors),NormFactors=norm_factors)
# write to a file that is useable
write.table(norm_factor_dt,norm_factor_file,sep="\t",quote=F,row.name=F)

