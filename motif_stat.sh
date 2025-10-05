#!/bin/bash
f_input=$1 #gid/space/motif
dis=$2 ; dis=${dis:=3} # default dis=3
dir=/scratch/users/yfyuan/PT/PTseq/D05_072922M4/PTgenomes32/txtFR

rm ${dir}/*txt_tmp

awk '{print $1,$2}' $f_input | while read gid motif ; do
  f_F=${dir}/${gid}_F.txt.${motif}.combine
  f_R=${dir}/${gid}_R.txt.${motif}.combine
  fout=${dir}/${gid}_pileup.txt

  # concat F, R txt files.
  awk 'FNR==NR{print $0" F"} FNR<NR{print $0" R"}' $f_F $f_R | sort -t ' ' -Vk1,2 >> ${fout}_tmp
done

echo -e "gid\ttotal\tC*AG\tC*CA\tC*CTG\tG*ATC\tG*AGC\tG*CTC\tG*AAC\tG*TTC"

awk '{print $1}' $f_input | sort -u | while read gid ; do
  fout=${dir}/${gid}_pileup.txt
  # input columns: 1scaffold, 2pos, 3cov, 4depR2, 5depR12, 6seq, 7strand
  echo -e "$gid\t\
  $(awk '($4>9.5||$5>9.5) {print $2}' ${fout}_tmp| wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......CAG/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......CCA/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......CCTG/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......GATC/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......GAGC/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......GCTC/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......GAAC/ {print $2}' ${fout}_tmp | wc -l)\t\
  $(awk '($4>9.5||$5>9.5)&&$6~/^......GTTC/ {print $2}' ${fout}_tmp | wc -l)"
done
  
