### in R version above 4
library("Rsubread")
options(stringsAsFactors = F)

input_dir = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/02-1_rmdup"
output_file = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/07_downstream_analysis/07.1_DEG/TTseq_featurecounts_count.txt"
output_log = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/07_downstream_analysis/07.1_DEG/TTseq_featurecounts.log"

annot_saf = read.table("/share/home/Blueberry/reference/annotation/ucsc/DLD1_cellline_INTAC/promoter/hg19.ucsc.refseq_proteincoding_fwd_and_rev_sense_INTAC.saf",
                      header = T,sep = "\t")


number_threads = 26
samplename_replace="_hg19.rmdup.q10.bam"


#annot_file="/share/home/Blueberry/reference/annotation/ucsc/DLD1_cellline_INTAC/hg19.ncbiRefSeq_DLD1.gtf"

files_use = dir(input_dir,pattern = paste0(samplename_replace,"$"),full.names = T)

# input_bam=c(file.path(input_dir, "TT-seq_dTAG-6h_MG132-3h_INTS11-8-dTAG_DLD1_rep2_211224_hg19.rmdup.q10.bam"),
#             file.path(input_dir, "TT-seq_dTAG-6h_MG132-3h_INTS11-8-dTAG_DLD1_rep1_211224_hg19.rmdup.q10.bam"))


### unstranded
# df_featurecount = featureCounts(files=input_bam,
#                                 annot.ext=annot_file,
#                                 isGTFAnnotationFile=TRUE,
#                                 GTF.featureType="transcript",
#                                 GTF.attrType="gene_id",
#                                 isPairedEnd=TRUE,
#                                 requireBothEndsMapped=TRUE,
#                                 countChimericFragments=FALSE,
#                                 countMultiMappingReads=FALSE)
# 
# ### reversely stranded
# df_featurecount_2 = featureCounts(files=input_bam,
#                                 annot.ext=annot_file,
#                                 isGTFAnnotationFile=TRUE,
#                                 GTF.featureType="transcript",
#                                 GTF.attrType="gene_id",
#                                 isPairedEnd=TRUE,
#                                 requireBothEndsMapped=TRUE,
#                                 countChimericFragments=FALSE,
#                                 countMultiMappingReads=FALSE,
#                                 strandSpecific=2)


df_featurecount_saf = featureCounts(files=files_use,
                                  annot.ext=annot_saf,
                                  isGTFAnnotationFile=FALSE,
                                  isPairedEnd=TRUE,
                                  requireBothEndsMapped=TRUE,
                                  countChimericFragments=FALSE,
                                  countMultiMappingReads=FALSE,
                                  strandSpecific=2,
                                  nthreads = number_threads)

df_result = df_featurecount_saf$counts
colnames(df_result) = gsub(samplename_replace,"",colnames(df_result) )

write.table(df_result, output_file,
      sep = "\t", quote = F, row.names = T, col.names = NA)

write.table(df_featurecount_saf$stat, output_log,
            sep = "\t", quote = F, row.names = T, col.names = NA)


