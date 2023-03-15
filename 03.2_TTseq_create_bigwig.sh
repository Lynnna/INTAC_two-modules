#!/usr/bin/bash



#================== bamcoverage(bw) ==================#
p_logs="$path_dir/logs"
p_rawdata="$path_dir/00-0_rawdata"
p_clean="$path_dir/01-0_cleandata"
p_cleanLog="${path_dir}/logs/01-0_cleandata"
p_align="$path_dir/02-0_align"
p_alignLog="$path_dir/logs/02-0_align"
p_rmdup="$path_dir/02-1_rmdup"
p_rmdupLog="$path_dir/logs/02-1_rmdup"
p_spike="$path_dir/02-2_spikein"
echo -e "\n
***************************
track of R1andR2 begins at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)
***************************"
p_track=$path_dir/03-0_track
p_trackLog=$path_dir/logs/03-0_track
mkdir -p ${p_track} ${p_trackLog}

cat ${p_spike}/scalefactor.txt | while read id;
do
arr=($id)
sample=${arr[0]}
scalefactor=${arr[10]}
bam_file=${p_rmdup}/${sample}_hg19.rmdup.Q7.bam
echo "=== sample: $sample ==="
    if [ ! -s "${p_track}/${sample}_rmdup-Q7_dr_fwd.bw" ]
    then
        echo "Generating file: ${p_track}/${sample}_rmdup-Q7_dr_fwd.bw..."
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${p_track}/${sample}_rmdup-Q7_dr_fwd.bw \
        --binSize 1 \
        --scaleFactor  $scalefactor \
        --numberOfProcessors 25 \
        --normalizeUsing None \
        --filterRNAstrand forward 
    else
        echo "${p_track}/${sample}_rmdup-Q7_dr_fwd.bw exists, continue..."
    fi
    
    if [ ! -s "${p_track}/${sample}_rmdup-Q7_dr_rev_minus.bw" ]
    then
    echo "Generating file: ${p_track}/${sample}_rmdup-Q7_dr_rev_minus.bw..."
        bamCoverage \
        --bam ${bam_file} \
        --skipNonCoveredRegions \
        --outFileName ${p_track}/${sample}_rmdup-Q7_dr_rev_minus.bw \
        --binSize 1 \
        --scaleFactor  -"${scalefactor}" \
        --numberOfProcessors 25 \
        --normalizeUsing None \
        --filterRNAstrand reverse
    else
    echo "${p_track}/${sample}_rmdup-Q7_dr_rev_minus.bw exists, continue"
    fi
done
