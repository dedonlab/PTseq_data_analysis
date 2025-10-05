#!/bin/bash
# set argvs
genome=$1
r1=$2
r2=$3
sample=$4 # job name

bow=$(basename ${genome} | cut -d '.' -f1) # name of bowtie db
bam="$sample".bam
bam_clean="$sample".clean.bam
t=1 # thread = 1

# dependencies
# bowtie2
# samtools
# bedtools
# SMARTcleaner-maste

# index reference genome
bowtie2-build -q --threads $t --seed 1 $genome $bow

samtools faidx $genome

# map
bowtie2 --quiet --sensitive --threads $t --seed 1 -x $bow -1 $r1 -2 $r2 | samtools sort -@ $t -O BAM -o $bam

# clean TTTTT off-target reads
SMARTcleaner cleanPEbam $genome $bam -o .

rm $sample.bam.noise.bam || true

bam_F=$sample'_clean_F.bam'
bam_R=$sample'_clean_R.bam'

samtools view -b -@ $t -f 163 -o ${bam_F} $bam_clean
samtools view -b -@ $t -f 147 -o ${bam_R} $bam_clean

dir_o=.

# identify pileups for R and F strand separately.
PAIRS=(F R)

for pair in ${PAIRS[*]} ; do
  bam_io=$dir_o/$sample'_clean_'$pair'.bam'
  cov_all_io=$dir_o/$sample'_all.'$pair'cov'
  cov_5_io=$dir_o/$sample'_5.'$pair'cov'
  cov_all_io_s=$dir_o/$sample'_all.'$pair'cov_sort'
  cov_5_io_s=$dir_o/$sample'_5.'$pair'cov_sort'
  cov_cmb=$dir_o/$sample'_cmb.'$pair'cov'
  # calculate coverage at n position across genome.
  bedtools genomecov -ibam $bam_io -d > $cov_all_io
  # calculate # of reads start at n position across genome.
  bedtools genomecov -ibam $bam_io -d -5 > $cov_5_io

  # combine 2 coverage files cover_all_io, cover_start_io.
  # 4 columns: 'scaffold','position','coverage_io','#read_start_here_io'
  join -j 2 -o 1.1,1.2,1.3,2.3 $cov_all_io $cov_5_io | awk -F' ' '{if ($3==0) print $0,0; else print $0,$4/$3;}' > $cov_cmb

  rm $cov_all_no $cov_all_io $cov_5_no $cov_5_io $cov_all_no_s $cov_all_io_s $cov_5_no_s $cov_5_io_s || true
done

# 5 columns:
# scaffold,  pos,  cov,  depth,  dep/cov_ratio

dir_w=.
Fcov=${dir_w}/${sample}_cmb.Fcov
Rcov=${dir_w}/${sample}_cmb.Rcov
Fpos=${dir_w}/${sample}_pileup_dep0_F.pos
Rpos=${dir_w}/${sample}_pileup_dep0_R.pos

awk '$4>0{print $0}' $Fcov > $Fpos
awk '$4>0{print $0}' $Rcov > $Rpos

# input: _pileup_dep0_F.pos _pileup_dep0_R.pos
# input columns: scaffold, pos, cov_at_pos, pileup_depth, depth/cov_ratio, sequence

# output: _pileup_dep0_F.pos.txt _pileup_dep0_R.pos.txt
# columns: scaffold, pos, cov_at_pos, pileup_depth, depth/cov_ratio, sequence (6 flank nt)

# retrieve 6nt flank sequences
f=6
# F.
seq_out=${Fpos}.txt
>${seq_out}.tmp

awk '{print $1,$2}' $Fpos | while read a b; do
  posl=$(($b-$f))
  if [ $posl -lt 0 ]; then
    posl=1
  fi
  posr=$(($b+$f))
  samtools faidx -c -i $genome ${a}:${posl}-${posr} | grep -v '^>' >> ${seq_out}.tmp
done

paste -d ' ' $Fpos ${seq_out}.tmp > ${seq_out}


# R.

seq_out=${Rpos}.txt
>${seq_out}.tmp
awk '{print $1,$2}' $Rpos | while read a b; do
  posl=$(($b-$f))
  if [ $posl -lt 0 ]; then
    posl=1
  fi
  posr=$(($b+$f))
  samtools faidx -c $genome ${a}:${posl}-${posr} | grep -v '^>' >> ${seq_out}.tmp
done
paste -d ' ' $plup ${seq_out}.tmp > ${seq_out}

# Merge file of F R strand.
# extract sequences from txt file to .fasta format.
awk '{print ">"$1"_"$2"\n"$6}' ${Fpos}.txt > ${sample}_pileup_dep0.fasta

awk '{print ">"$1"_"$2"\n"$6}' ${Rpos}.txt >> ${sample}_pileup_dep0.fasta

