options(stringsAsFactors = F)
library("DESeq2")
library("tidyverse")
library(edgeR)
library(ggplot2)

input_file = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/07_downstream_analysis/07.1_DEG/TTseq_featurecounts_count.txt"
spike_in = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/02-2_spikein/scalefactor.txt"
output_path = "/share/home/Blueberry/Projects/Chenlab/Hushibin/INTAC/TT-seq/07_downstream_analysis/07.1_DEG"
output_file = file.path(output_path,"Deseq2_result")

# count matrix
df_matrix = read.table(input_file,header = T,sep = "\t",row.names = 1)
colnames(df_matrix) = gsub("DLD1_","",gsub("TT.seq_","",gsub("_21....$","",colnames(df_matrix))))

#scale factor 
df_spike = read.table(spike_in,header = F, sep = "\t")
row.names(df_spike) = df_spike$V1
row.names(df_spike) = gsub("-",".",gsub("DLD1_","",gsub("TT-seq_","",gsub("_21....$","",row.names(df_spike) ))))

# metadata
control_1 = grep("rep1",grep("^DMSO",sort(colnames(df_matrix)),value = T), value = T)
control_2 = grep("rep2",grep("^DMSO",sort(colnames(df_matrix)),value = T), value = T)

treatment_1 = grep("rep1",grep("^dTAG",sort(colnames(df_matrix)),value = T), value = T)
treatment_2 = grep("rep2",grep("^dTAG",sort(colnames(df_matrix)),value = T), value = T)

df_meta=data.frame(control_1 = control_1,
                   control_2 = control_2,
                   treatment_1 = treatment_1,
                   treatment_2 = treatment_2)

#i = 1
df_combine = NULL
for ( i in 1:nrow(df_meta)){

meta_use = df_meta[i,] %>% 
  gather(key = 'condition', value = 'sample')

row.names(meta_use) = meta_use$sample
meta_use$group = sapply(strsplit(meta_use$condition,"_"),"[",1)
meta_use$group = factor(meta_use$group,levels = c("control","treatment") )
meta_use = meta_use[,!grepl("sample",colnames(meta_use))]
scale_use = df_spike[row.names(meta_use),11]
names(scale_use) = row.names(meta_use)

# QC filter genes
matrix_use = df_matrix[,row.names(meta_use)]
matrix_use = sweep(matrix_use, MARGIN = 2, scale_use, `*`) 
matrix_use = round(matrix_use)
matrix_use = matrix_use[rowSums(cpm(matrix_use) >1) >=2,]
  

# build directory and get sample name
samplename <- gsub("_rep1","", paste0(row.names(meta_use)[meta_use$condition=="treatment_1"],
                     "_vs_",
                     row.names(meta_use)[meta_use$condition=="control_1"]))
output_path_use <- file.path(output_file,samplename)

if (!file.exists(output_path_use)){
  dir.create(output_path_use)
  }

write.table(matrix_use,
            file.path(output_path_use,paste0(samplename,'_adjusted_count_matrix.txt')), sep='\t', quote=FALSE,col.names = NA,row.names = T)



# The sizeFactors vector assigns to each column of the count matrix a value, the size factor, 
# such that count values in the columns can be brought to a common scale by dividing by the 
# corresponding size factor (as performed by counts(dds, normalized=TRUE)).
# 我们之前计算的是乘，现在是除，所以要1/raw_scale_factor
# 
# scale_use = 1/df_spike[row.names(meta_use),11]
# names(scale_use) = row.names(meta_use)


if (all(rownames(meta_use) == colnames(matrix_use)) == TRUE){
  dds <- DESeqDataSetFromMatrix(countData = matrix_use,
                                colData = meta_use,
                                design = ~ group)
  # sizeFactors(dds) <- scale_use
  sizeFactors(dds) <- c(1,1,1,1)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("group","treatment","control"))
  baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$group == "treatment"])
  baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$group == "control"])
  # unique(names(baseMeanA)==row.names(res))
  res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
  resOrdered <- res[order(res$pvalue),]
  
  write.table(as.data.frame(resOrdered),
              file.path(output_path_use,paste0(samplename,'_DESeq2.DE_results.txt')), sep='\t', quote=FALSE,col.names = NA,row.names = T)
  
  # diff_gene_deseq2 <-as.data.frame(subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
  
  
  #DEG

  DEG_num <- nrow(resOrdered[resOrdered$padj <0.05,])
  log2FC_up <- nrow(subset(resOrdered,padj < 0.05 &  log2FoldChange > 1))
  log2FC_down <- nrow(subset(resOrdered,padj < 0.05 &  log2FoldChange < -1))
  sample_treat <-  gsub("_rep1","",row.names(meta_use)[meta_use$condition=="treatment_1"])
  sample_con <-    gsub("_rep1","",row.names(meta_use)[meta_use$condition=="control_1"])
  
  com <-c(sample_treat,sample_con,DEG_num,log2FC_down,log2FC_up)
  df_combine <- rbind(df_combine,com)
  
  
  
  ## data for plotting
  res_DMSO5 <- res[res$baseMeanB>= 0,]
  res_DMSO5$color <- ifelse(res_DMSO5$padj<0.05, "Differential", "Not-Differential")
  res_DMSO5$color[res_DMSO5$color == "Differential"] <- ifelse(res_DMSO5[res_DMSO5$color == "Differential",]$log2FoldChange > 0,
                                                               "Up-Regulated", "Down-Regulated")
  
  #### MA plot ####
  title <- sprintf("Number of features = %s\nUp-Regulated = %s (%s Percent)\nDown-Regulated = %s (%s Percent)",
                   nrow(res_DMSO5), sum(res_DMSO5$color=="Up-Regulated"), 
                   round(100 * sum(res_DMSO5$color=="Up-Regulated") / nrow(res_DMSO5), 2),
                   sum(res_DMSO5$color=="Down-Regulated"), 
                   round(100 * sum(res_DMSO5$color=="Down-Regulated") / nrow(res_DMSO5), 2)
  )
  
  pal <- c("Up-Regulated" = "firebrick3", "Not-Differential" = "lightgrey", "Down-Regulated" = "dodgerblue3")
  

  plot1=ggplot(res_DMSO5,
               aes(log2(baseMean+1),log2FoldChange))+
    geom_point(aes(colour=color), size = 0.7 ) + scale_colour_manual(values = pal) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme_classic() +
    labs(x = "Log2 Mean",y="Log2 Fold Change",title = title ) +
    theme(axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          axis.ticks =element_line(colour = "black"))+
    theme(aspect.ratio=1)+
    guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(filename = file.path(output_path_use,paste0(samplename,'_DESeq2_MAplot.pdf')),
         plot1, width = 5,height = 5)
  
 ####Volcano Plot####
  plot2= ggplot(res_DMSO5,
                aes(log2FoldChange,-1*log10(padj)))+
    geom_point( aes(colour=color),size = 0.7) + 
    scale_colour_manual(values = pal) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme_classic() +
    labs(x = "log2FC",y="-1*log10(FDR)",title = title ) +
    theme(axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          axis.ticks =element_line(colour = "black"))+
    theme(aspect.ratio=1)+
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  ggsave(filename = file.path(output_path_use,paste0(samplename,'_DESeq2_VolcanoPlot.pdf')),
         plot2, width = 5,height = 5)
  
  plot3= ggplot(res_DMSO5,
                aes(log2FoldChange,-1*log10(padj)))+
    geom_point( aes(colour=color),size = 0.7) + 
    scale_colour_manual(values = pal) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme_classic() +
    labs(x = "log2FC",y="-1*log10(FDR)",title = title ) +
    theme(axis.text.x = element_text(colour="black"), 
          axis.text.y = element_text(colour="black"),
          axis.ticks =element_line(colour = "black"))+
    theme(aspect.ratio=1)+
    ylim(c(0,50))+xlim(c(-4,4))+
    guides(color = guide_legend(override.aes = list(size = 3)))
    
  ggsave(filename = file.path(output_path_use,paste0(samplename,'_DESeq2_VolcanoPlot_changebar.pdf')),
         plot3, width = 5,height = 5)
  }

}
  
colnames(df_combine) <- c("sample_treat","sample_con","DEG_num","log2FC_down(<-1)","log2FC_up(>1)")

write.table(df_combine, file.path(output_path,'TTseq_DESeq2.DE_results_total.txt'), 
            sep='\t', quote=FALSE,col.names = T,row.names = F)
  
  
