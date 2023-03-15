####################################################################################################
#########################   ChIP-seq analysis pipeline  with spike-in     ##########################
####################################################################################################

###load modules
module load genomics/apps
module load miniconda/3
source activate chipseq



###file_path
input_path=/ChIP-seq/data_unchanged
file_path=/ChIP-seq
raw_dir=${file_path}/00_rawdata
logs_dir=${file_path}/logs
fastqc_dir=${file_path}/01_rawfastqc

trimmedFastq_dir=${file_path}/02_trimmeddata
trimmedFastq_log_dir=${file_path}/logs/trimmeddata

trimmeddata_fastqc_dir=${file_path}/02_trimmeddata_fastqc


align_exp_dir=${file_path}/03_bam_hg19
alignexp_log_dir=${file_path}/logs/align_hg19

align_spike_dir=${file_path}/03_spikebam_mm10
alignspike_log_dir=${file_path}/logs/align_mm10

exp_bam_rmdup=${file_path}/03_bam_hg19_rmdup
rmdup_exp_log=${file_path}/logs/rmdup_state_hg19

spike_bam_rmdup=${file_path}/03_spikebam_mm10_rmdup
rmdup_spike_log=${file_path}/logs/rmdup_state_mm10


bw_fulllength_dir=${file_path}/04_bw_fulllength

peak_dir=${file_path}/05_peak
peak_log_dir=${file_path}/logs/callpeak

peak_alone_dir=${file_path}/05_peak/05.1_peak_alone
peak_compare_dir=${file_path}/05_peak/05.2_peak_compare

plot_dir=${file_path}/06_analysis
sampleinfo=${file_path}/sampleinfo.txt

sampledir=${file_path}/Sampleinfo

#reference genome
GENOME_EXP="/share/home/Blueberry/reference/index/bowtie2/hg19/hg19"
GENOME_SPIKE="/share/home/Blueberry/reference/index/bowtie2/bowtie2_download/mm10/mm10"
blacklist="/share/home/Blueberry/reference/annotation/encode/blacklist/hg19-blacklist.v2.bed"


###create  sampleinfo

mkdir -p ${sampledir}
cat $sampleinfo |cut -f 2 > ${sampledir}/sample_use.txt

cat ${sampledir}/sample_use.txt|grep -i "input"|while read id
do
    input_sample=${id/*Input_/} 
    echo $id > ${sampledir}/IP_${input_sample}.txt
    grep ${input_sample} ${sampledir}/sample_use.txt | grep -v -i "input" >> ${sampledir}/IP_${input_sample}.txt
done

head -n 1 ${sampledir}/IP*txt | grep "^IP" > ${sampledir}/input.txt


#step 0
####check files 
mkdir -p $logs_dir
for id in `ls ${input_path}/*/md5.txt`
do
path1=$(dirname $id)
echo $path1
cd $path1
md5sum -c $id
done > ${logs_dir}/checkfiles.log


#step 1.1
####change file name#####

mkdir -p ${raw_dir}

cat $sampleinfo| while read id;
do
arr=($id)
sample1=${arr[0]}
sample2=${arr[1]}
fq1=$(ls ${input_path}/Cleandata/*/*R1.fq.gz|grep "${sample1}_R1.fq.gz")
fq2=$(ls ${input_path}/Cleandata/*/*R2.fq.gz|grep "${sample1}_R2.fq.gz")
ln -s $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
ln -s $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
done

#step 1.2
####fastqc of raw data ####
mkdir -p ${fastqc_dir}

in_path=${raw_dir}
out_path=${fastqc_dir}

nohup_number=0
for file in `ls ${in_path}/*R1.fastq.gz`
do
    ID=$(basename $file)
    fq1=${file}
    fq2=${file/R1/R2}
    if [ ! -s ${out_path}/"${ID/.fastq.gz/_fastqc.zip}" ]
    then
    echo "Generating file: $path/00-1_rawFastqc/"${ID}_R1_fastqc.zip"..."
    fastqc $fq1 -t 1  -o ${out_path}/  &
    fastqc $fq2 -t 1  -o ${out_path}/  &
    nohup_number=`echo $nohup_number+2 | bc`
    fi
    
    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done

wait

#step 1.3
### merge reports of fastqc
multiqc ${fastqc_dir}/ -n rawdata_multiqc -o ${fastqc_dir}/


#step 2.1
### Trimming adapters  (trim_galore)

mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}

nohup_number=0
for fq1 in `ls ${raw_dir}/*R1.fastq.gz`
do
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
    if [ ! -s ${trimmedFastq_dir}/"$(basename ${fq1/R1.fastq.gz/trimmed_R1.fq.gz})" ]
    then
     trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired -o ${trimmedFastq_dir} $fq1 $fq2 \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fq.gz/_trimmed.log})" 2>&1 &

     nohup_number=`echo $nohup_number+1 | bc`
    fi

   
    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done

wait

#step 2.2
### Cleaning up filenames in trimmedFastq (bowtie automatically names PE --un output)

for FILE in ${trimmedFastq_dir}/*1.fq.gz
do
    if [ ! -s ${FILE/R1_val_1.fq.gz/trimmed_R1.fq.gz} ]
    then
        mv "$FILE" ${FILE/R1_val_1.fq.gz/trimmed_R1.fq.gz}
    fi
done

for FILE in ${trimmedFastq_dir}/*2.fq.gz
do 
    if [ ! -s ${FILE/R2_val_2.fq.gz/trimmed_R2.fq.gz} ]
    then
        mv "$FILE" ${FILE/R2_val_2.fq.gz/trimmed_R2.fq.gz}
    fi
done

#step 2.3
#####qc for trimmed data####
mkdir -p ${trimmeddata_fastqc_dir}


nohup_number=0
for file in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do
    ID=$(basename $file)
    fq1=${file}
    fq2=${file/R1/R2}
    if [ ! -s ${trimmeddata_fastqc_dir}/"${ID/.fq.gz/_fastqc.html}" ]
    then
    echo "Generating file: ${trimmeddata_fastqc_dir}/"${ID/.fq.gz/_fastqc.html}"..."
    fastqc $fq1 -t 1  -o ${trimmeddata_fastqc_dir}/  &
    fastqc $fq2 -t 1  -o ${trimmeddata_fastqc_dir}/  &

     nohup_number=`echo $nohup_number+2 | bc`
    fi
   
    if [[ $nohup_number -eq 28 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
done



wait

### merge reports of fastqc
multiqc ${trimmeddata_fastqc_dir}/ -n trimmeddata_multiqc -o ${trimmeddata_fastqc_dir}/



#step 3.1
###Aligning to experimental genome#####
echo -e "\n***************************\nalign of experimental genome begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${align_exp_dir}
mkdir -p ${alignexp_log_dir}
mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}


for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do 
fq2=${fq1/R1.fq.gz/R2.fq.gz}
sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
if [ ! -s ${exp_bam_rmdup}/${sample}_hg19.rmdup.bam ]
then
    (bowtie2  -p 25   -x  ${GENOME_EXP} -N 1  -1 $fq1 -2 $fq2 \
    2> ${alignexp_log_dir}/${sample}_align.log) \
    |samtools view -bS -F 3844 -f 2 -q 30 \
    |samtools sort  -O bam  -@ 25 -o ${align_exp_dir}/${sample}_hg19.bam 
    samtools index ${align_exp_dir}/${sample}_hg19.bam
     
    # picard

    picard MarkDuplicates -REMOVE_DUPLICATES True \
        -I ${align_exp_dir}/${sample}_hg19.bam \
        -O ${exp_bam_rmdup}/${sample}_hg19.rmdup.bam \
        -M ${exp_bam_rmdup}/${sample}_hg19.rmdup.metrics
    samtools index  ${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
    samtools flagstat ${exp_bam_rmdup}/${sample}_hg19.rmdup.bam > ${rmdup_exp_log}/${sample}_hg19.rmdup.stat
    
fi
done 


#step 3.2
###Aligning to spike-in genome#####
echo -e "\n***************************\nalign of spike-in genome begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"

mkdir -p ${align_spike_dir}
mkdir -p ${alignspike_log_dir}
mkdir -p ${spike_bam_rmdup}
mkdir -p ${rmdup_spike_log}


for fq1 in `ls ${trimmedFastq_dir}/*R1.fq.gz`
do 
fq2=${fq1/R1.fq.gz/R2.fq.gz}
sample="$(basename ${fq1/_trimmed_R1.fq.gz/})"
if [ ! -s ${spike_bam_rmdup}/${sample}_onlymm10.rmdup.bam ]
then
    (bowtie2  -p 25   -x  ${GENOME_SPIKE} -N 1  -1 $fq1 -2 $fq2 \
    2> ${alignspike_log_dir}/${sample}_spikealign.log) \
    |samtools view -bS -F 3844 -f 2 -q 30 \
    |samtools sort  -O bam  -@ 25 -o ${align_spike_dir}/${sample}_onlymm10.bam 
    samtools index ${align_spike_dir}/${sample}_onlymm10.bam
     
    # picard

    picard MarkDuplicates -REMOVE_DUPLICATES True \
        -I ${align_spike_dir}/${sample}_onlymm10.bam \
        -O ${spike_bam_rmdup}/${sample}_onlymm10.rmdup.bam \
        -M ${spike_bam_rmdup}/${sample}_onlymm10.rmdup.metrics
    samtools index  ${spike_bam_rmdup}/${sample}_onlymm10.rmdup.bam
    samtools flagstat ${spike_bam_rmdup}/${sample}_onlymm10.rmdup.bam > ${rmdup_spike_log}/${sample}_onlymm10.rmdup.stat
    
fi
done 





#step 3.3
### calculate normalization factors ###
echo -e "\n***************************\nCalculating normalization factors at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n***************************"


hg19_path=${rmdup_exp_log}
mm10_path=${rmdup_spike_log}
align_path=${alignexp_log_dir}

if [ -s "${logs_dir}/scalefactor.txt" ] 
then
    rm  ${logs_dir}/scalefactor.txt
fi

touch ${logs_dir}/scalefactor.txt

if [ -s "${logs_dir}/spike-input.txt" ] 
then
    rm  ${logs_dir}/spike-input.txt
fi

touch ${logs_dir}/spike-input.txt

# input file
input_file=${sampledir}/input.txt
spike_sample=`head -n 1 ${input_file}`

spike_hgRatio=$(cat $hg19_path/$spike_sample"_hg19.rmdup.stat" | grep "total (QC-passed reads"|cut -d " " -f 1)
spike_mmRatio=$(cat $mm10_path/$spike_sample"_onlymm10.rmdup.stat" | grep "total (QC-passed reads"|cut -d " " -f 1)
spike_product=$(echo $spike_hgRatio/$spike_mmRatio | bc -l)

cat $input_file | while read sample;do
    total_log=${align_path}/${sample}_align.log
    mm10_log=${mm10_path}/${sample}_onlymm10.rmdup.stat
    hg19_log=${hg19_path}/${sample}_hg19.rmdup.stat

    ALLREADS=$(cat ${total_log}|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    Hg19_READS=$(cat ${total_log}| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 | \
    sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    Hg19_RATIO=$(cat ${total_log}|grep "overall alignment rate"|cut -d "%" -f 1)

    Hg19_qc_READS=$(cat ${hg19_log}|grep "total (QC-passed reads"|cut -d " " -f 1)
    Hg19_qc_RATIO=$(echo "${Hg19_qc_READS}/${ALLREADS}"|bc -l)

    MM10_qc_READS=$(cat ${mm10_log}|grep "total (QC-passed reads"|cut -d " " -f 1)
    QC_reads=$(echo "${MM10_qc_READS}+${Hg19_qc_READS}"|bc )
    MM10_qc_RATIO_intotal=$(echo "${MM10_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_RATIO_inqc=$(echo "${MM10_qc_READS}/${QC_reads}"|bc -l)
    spike_factor=$(echo "($Hg19_qc_READS/$MM10_qc_READS)/$spike_product"|bc -l)

    echo -e $sample"\t"$ALLREADS"\t"$Hg19_READS"\t"$Hg19_RATIO"\t"$Hg19_qc_READS"\t"$Hg19_qc_RATIO"\t"$MM10_qc_READS"\t"$QC_reads"\t"$MM10_qc_RATIO_intotal"\t"$MM10_qc_RATIO_inqc"\t"$spike_factor  >> ${logs_dir}/spike-input.txt
done


for file in ${sampledir}/IP*txt;
do
    input=`head -n 1 $file`
    spike_factor=$(grep $input ${logs_dir}/spike-input.txt | cut -f 11)
    cat $file|while read sample;do

    total_log=${align_path}/${sample}_align.log
    mm10_log=${mm10_path}/${sample}_onlymm10.rmdup.stat
    hg19_log=${hg19_path}/${sample}_hg19.rmdup.stat

    ALLREADS=$(cat ${total_log}|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
    Hg19_READS=$(cat ${total_log}| sed 's/%//g' | awk '{printf $0"\t"}'  |cut -f 4,5,8,13,14 | \
    sed 's/\t/\n/g' | awk '{print $1}' | awk '{printf $0"\t"}'|awk '{print 2*($1+$2+$3)+$4+$5}')
    Hg19_RATIO=$(cat ${total_log}|grep "overall alignment rate"|cut -d "%" -f 1)

    Hg19_qc_READS=$(cat ${hg19_log}|grep "total (QC-passed reads"|cut -d " " -f 1)
    Hg19_qc_RATIO=$(echo "${Hg19_qc_READS}/${ALLREADS}"|bc -l)

    MM10_qc_READS=$(cat ${mm10_log}|grep "total (QC-passed reads"|cut -d " " -f 1)
    QC_reads=$(echo "${MM10_qc_READS}+${Hg19_qc_READS}"|bc )
    MM10_qc_RATIO_intotal=$(echo "${MM10_qc_READS}/${ALLREADS}"|bc -l)
    MM10_qc_RATIO_inqc=$(echo "${MM10_qc_READS}/${QC_reads}"|bc -l)
    SCALEFACTOR=$(echo "1/((${MM10_qc_READS}/1000000)*$spike_factor)" | bc -l )

    echo -e $sample"\t"$ALLREADS"\t"$Hg19_READS"\t"$Hg19_RATIO"\t"$Hg19_qc_READS"\t"$Hg19_qc_RATIO"\t"$MM10_qc_READS"\t"$QC_reads"\t"$MM10_qc_RATIO_intotal"\t"$MM10_qc_RATIO_inqc"\t"$SCALEFACTOR >> ${logs_dir}/scalefactor.txt
    done
done



wait


#step 4.1
### Making RPKM-normalized bigWig files with full-length reads without spike-in ###
echo -e "\n***************************\ntracks need to be done....\n***************************"

mkdir -p ${bw_fulllength_dir}

cat  ${logs_dir}/scalefactor.txt | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[10]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam

    if [ ! -s "${bw_fulllength_dir}/${sample}_fulllength.bw" ]
    then
        bamCoverage -b ${bam_file} \
        --binSize 1 \
        --blackListFileName  ${blacklist} \
        --normalizeUsing None \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 23 \
        -o ${bw_fulllength_dir}/${sample}_fulllength.bw 2>${bw_fulllength_dir}/${sample}.log
    fi
done


      
#step 5.1
###call peak
echo -e "\n***************************\ncall-peak begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S); control-IP-peak need to be done....\n***************************"

mkdir -p ${peak_log_dir}
mkdir -p ${peak_alone_dir}
mkdir -p ${peak_compare_dir}

####without input
nohup_number=0
cat  ${logs_dir}/scalefactor.txt | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[10]}
bam_file=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam
if [ ! -s ${peak_alone_dir}/${sample}_alone_peaks.narrowPeak ]
then
    macs2 callpeak -t $bam_file -f BAMPE -n ${sample}_alone -g hs --keep-dup all --outdir ${peak_alone_dir} \
    2> ${peak_log_dir}/${sample}_alone.log &

    nohup_number=`echo $nohup_number+1 | bc`
        fi

    
    if [[ $nohup_number -eq 2 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi

done

wait

nohup_number=0
cat  ${logs_dir}/scalefactor.txt | while read id;
do
arr=($id)
sample=${arr[0]}
if [ ! -s ${peak_alone_dir}/${sample}_alone_peaks.final.narrowPeak ]
then
    bedtools intersect -a ${peak_alone_dir}/${sample}_alone_peaks.narrowPeak \
    -b $blacklist \
    -f 0.25 -v > ${peak_alone_dir}/${sample}_alone_peaks.final.narrowPeak &


    nohup_number=`echo $nohup_number+1 | bc`
fi

    if [[ $nohup_number -eq 27 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi

done

wait



####with input as control
nohup_number=0
for file in ${sampledir}/IP*txt;
do
    input=`head -n 1 $file`
    input_bam=${exp_bam_rmdup}/${input}_hg19.rmdup.bam
    sed 1d $file|cat |while read sample ;
    do
        if [ ! -s ${peak_compare_dir}/${sample}_compare_peaks.narrowPeak ]
        then 
            IP_bam=${exp_bam_rmdup}/${sample}_hg19.rmdup.bam 
            macs2 callpeak -c $input_bam \
                           -t $IP_bam \
                           -f BAMPE -n ${sample}_compare -g hs --keep-dup all --outdir ${peak_compare_dir} \
                           2> ${peak_log_dir}/${sample}_compare.log &

            nohup_number=`echo $nohup_number+1 | bc`
        fi

    
    if [[ $nohup_number -eq 2 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi
  done
done

wait


nohup_number=0
cat  ${logs_dir}/scalefactor.txt | while read id;
do
arr=($id)
sample=${arr[0]}
if [ ! -s ${peak_compare_dir}/${sample}_compare_peaks.final.narrowPeak ]
then
    bedtools intersect -a ${peak_compare_dir}/${sample}_compare_peaks.narrowPeak \
    -b $blacklist \
    -f 0.25 -v > ${peak_compare_dir}/${sample}_compare_peaks.final.narrowPeak &
    nohup_number=`echo $nohup_number+1 | bc`
fi


    if [[ $nohup_number -eq 27 ]]
    then
        wait
        echo "waiting..."
        nohup_number=0
    fi

done

wait














