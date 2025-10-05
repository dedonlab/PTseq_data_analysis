

genome=2807EA_1118_063_H5_final.scaffolds.fasta
bow=2807EA_1118_063_H5_bowtie
sample=Lanchnosp
bam="$sample".bam
bam_clean="$sample".clean.bam
t=1

module load bowtie2/2.4.5
module load samtools/1.19.2
module load bedtools/2.30.0

bowtie2-build -q --threads $t --seed 1 $genome $bow

samtools faidx $genome

r1=220214Ded_D22-1576_1_final.fastq
r2=220214Ded_D22-1576_2_final.fastq

bowtie2 --quiet --sensitive --threads $t --seed 1 -x $bow -1 $r1 -2 $r2 | samtools sort -@ $t -O BAM -o $bam

export PATH=/blue/lagard/yuanyifeng/software/SMARTcleaner-master:$PATH
SMARTcleaner cleanPEbam $genome $bam -o .

rm $sample.bam.noise.bam || true

bam_F=$sample'_clean_F.bam'
bam_R=$sample'_clean_R.bam'

samtools view -b -@ $t -f 163 -o ${bam_F} $bam_clean
samtools view -b -@ $t -f 147 -o ${bam_R} $bam_clean

dir_o=.

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
  if [[ $(basename $genome) == *'_final.scaffolds.fasta' ]] ; then
    # for genomes in Alm's collection.
    awk '{split($1, subfield, "_"); print subfield[1]"\t"$0}' $cov_all_io | sed 's/^scaffold//g' | sort -k1,1n -k3,3n | awk '{print $2"\t"$3"\t"$4}' > $cov_all_io_s
    awk '{split($1, subfield, "_"); print subfield[1]"\t"$0}' $cov_5_io | sed 's/^scaffold//g' | sort -k1,1n -k3,3n | awk '{print $2"\t"$3"\t"$4}' > $cov_5_io_s
    # combine 2 coverage files cover_all_io, cover_start_io.
    # 4 columns: 'scaffold','position','coverage_io','#read_start_here_io'.
    join -j 2 -o 1.1,1.2,1.3,2.3 $cov_all_io_s $cov_5_io_s | awk -F' ' '{if ($3==0) print $0,0; else print $0,$4/$3;}' > $cov_cmb
  elif [[ $(basename $genome) == 'GUT_'*'.fasta' ]] ; then
    # for genomes in UHGG's collection.
    awk '{split($1, subfield, "_"); print subfield[3]"\t"$0}' $cov_all_io | sort -k1,1n -k3,3n | awk '{print $2"\t"$3"\t"$4}' > $cov_all_io_s
    awk '{split($1, subfield, "_"); print subfield[3]"\t"$0}' $cov_5_io | sort -k1,1n -k3,3n | awk '{print $2"\t"$3"\t"$4}' > $cov_5_io_s
    # combine 2 coverage files cover_all_io, cover_start_io.
    # 4 columns: 'scaffold','position','coverage_io','#read_start_here_io'
    join -j 2 -o 1.1,1.2,1.3,2.3 $cov_all_io_s $cov_5_io_s | awk -F' ' '{if ($3==0) print $0,0; else print $0,$4/$3;}' > $cov_cmb
  else
    #sort -k1,1n -k2,2n $cov_all_io | awk '{print $1"\t"$2"\t"$3}' > $cov_all_io_s
    #sort -k1,1n -k2,2n $cov_5_io | awk '{print $1"\t"$2"\t"$3}' > $cov_5_io_s
    join -j 2 -o 1.1,1.2,1.3,2.3 $cov_all_io $cov_5_io | awk -F' ' '{if ($3==0) print $0,0; else print $0,$4/$3;}' > $cov_cmb
  fi
  #rm $cov_all_no $cov_all_io $cov_5_no $cov_5_io $cov_all_no_s $cov_all_io_s $cov_5_no_s $cov_5_io_s || true
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
# input columns: scaffold, pos, cov_at_pos, pileup_dep_mapper_R2, pileup_dep_mapper_R12.
# scaffold,  pos,  cov,  depth,  dep/cov_ratio
#cutoff=0.5 #cutoff = 0.5.

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
