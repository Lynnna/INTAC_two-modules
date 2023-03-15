#!/usr/bin/bash


#================== align ==================#
p_logs="$path_dir/logs"
p_rawdata="$path_dir/00-0_rawdata"
p_clean="$path_dir/01-0_cleandata"
p_cleanLog="${path_dir}/logs/01-0_cleandata"
echo -e "\n
***************************
Align begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
p_align="$path_dir/02-0_align"
p_alignLog="$path_dir/logs/02-0_align"
mkdir -p ${p_align} ${p_alignLog}

index=/share/home/Blueberry/reference/index/star/hg19_star_index
mm10index=/share/home/Blueberry/reference/index/star/mm10_star_index
cat $path_dir/samples | while read i;do
    CheckBam ${p_align}/${i}_hg19.Aligned.sortedByCoord.out.bam
    if [ ! -s ${p_align}/${i}_hg19.Aligned.sortedByCoord.out.bam ]
    then
        STAR --runThreadN 27 --genomeDir $index \
        --readFilesCommand zcat --readFilesIn ${p_clean}/${i}_R1_val_1.fq.gz ${p_clean}/${i}_R2_val_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix ${p_align}/${i}_hg19.
    fi
done
cat $path_dir/samples | while read i;do
    CheckBam ${p_align}/${i}_mm10.Aligned.sortedByCoord.out.bam
    if [ ! -s ${p_align}/${i}_mm10.Aligned.sortedByCoord.out.bam ]
    then
    STAR --runThreadN 27 --genomeDir $mm10index \
        --readFilesCommand zcat --readFilesIn ${p_clean}/${i}_R1_val_1.fq.gz ${p_clean}/${i}_R2_val_2.fq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --outFilterType BySJout --outFilterMultimapNmax 1 \
        --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 \
        --alignIntronMin 20 --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix ${p_align}/${i}_mm10.
    fi
done


#================== rmdup ==================#
p_logs="$path_dir/logs"
p_rawdata="$path_dir/00-0_rawdata"
p_clean="$path_dir/01-0_cleandata"
p_cleanLog="${path_dir}/logs/01-0_cleandata"
p_align="$path_dir/02-0_align"
p_alignLog="$path_dir/logs/02-0_align"
echo -e "\n
***************************
rmdup begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
p_rmdup="$path_dir/02-1_rmdup"
p_rmdupLog="$path_dir/logs/02-1_rmdup"
mkdir -p ${p_rmdup} ${p_rmdupLog}

cat $path_dir/samples | while read i;do
    if [ ! -s ${p_rmdupLog}/${i}_hg19.rmdup.Q7.stat ]
    then

    input=`ls ${p_align}/${i}_hg19*sortedByCoord.out.bam`
    rmdup=${p_rmdup}/${i}_hg19.rmdup.bam
    metrics=${p_rmdup}/${i}_hg19.metrics
    picard MarkDuplicates REMOVE_DUPLICATES=True \
        INPUT=$input OUTPUT=${rmdup} METRICS_FILE=${metrics} \
        2>${p_rmdupLog}/${i}_hg19.dup.log;
    rmdupQ7=${p_rmdup}/${i}_hg19.rmdup.Q7.bam
    samtools view -b -q 7 -f 2 ${rmdup} -o ${rmdupQ7}
    samtools index -@ 26 ${rmdupQ7}
    samtools flagstat -@ 25 ${rmdupQ7} > ${p_rmdupLog}/${i}_hg19.rmdup.Q7.stat
    echo "rmdupQ7 flagstat has done;file has generated in ${p_rmdupLog}/${i}_hg19.rmdup.Q7.stat"
    else
    echo "${p_rmdupLog}/${i}_hg19.rmdup.Q7.stat exists, continue..."
    fi

    if [ ! -s ${p_rmdupLog}/${i}_mm10.rmdup.Q7.stat ]
    then
    input=`ls ${p_align}/${i}_mm10*.sortedByCoord.out.bam`
    rmdup=${p_rmdup}/${i}_mm10.rmdup.bam
    metrics=${p_rmdup}/${i}_mm10.metrics
    picard MarkDuplicates REMOVE_DUPLICATES=True \
        INPUT=$input OUTPUT=${rmdup} \
        METRICS_FILE=${metrics} 2>${p_rmdupLog}/${i}_mm10.dup.log;
    samtools index -@ 25 ${rmdup}
    rmdupQ7=${p_rmdup}/${i}_mm10.rmdup.Q7.bam
    samtools view -b -q 7 -f 2 ${rmdup} -o ${rmdupQ7}
    samtools flagstat -@ 25 ${rmdupQ7} > ${p_rmdupLog}/${i}_mm10.rmdup.Q7.stat
    echo "rmdupQ7 flagstat has done;file has generated in ${p_rmdupLog}/${i}_mm10.rmdup.Q7.stat"
    else
    echo "${p_rmdupLog}/${i}_mm10.rmdup.Q7.stat exists, continue..."
    fi
done




