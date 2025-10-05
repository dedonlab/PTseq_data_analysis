#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanyifeng@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#############################################
module load bowtie2/2.4.5
module load samtools/1.19.2
genome=2807EA_1118_063_H5_final.scaffolds.fasta
bow=2807EA_1118_063_H5_bowtie
sample=Lanchnosp

dir_w=.
Fpos=${dir_w}/${sample}_pileup_dep0_F.pos
Rpos=${dir_w}/${sample}_pileup_dep0_R.pos

f=6
# F.
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

#R.
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
paste -d ' ' $Rpos ${seq_out}.tmp > ${seq_out}
