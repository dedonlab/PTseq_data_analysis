genome=2807EA_1118_063_H5_final.scaffolds.fasta
gff=2807EA_1118_063_H5_prokka.gff
sample=Lanchnosp
intergene_gff=${gff}.intergene

python3 gff_intergene.py ${gff} scaffold
awk -F '\t' '{print $3}' ${gff}.intergene | sort | uniq -c

motif=GAGC
f_raw="${sample}"_all_"$motif".txt
cat "$genome" | seqkit locate -i -p "$motif" > ${f_raw}.tmp

# re-format
# input:
# seqID patternName pattern strand start end matched
# scaffold1_size270648 GAAG GAAG + 174 177 GAAG
# scaffold1_size270648 GAAG GAAG + 433 436 GAAG

# desired output:
# 'scaffold','pileup_pos', 'cov','depth','ratio','sequence','strand','motif'
# scaffold1_size270648 1195 103 36 0.349515 TCTGTTGAAGGTA R GAAG
awk -F '\t' '{if ($4=="+") print $1, $5, "c d r s R", $2; else if($4=="-") print $1, $5, "c d r s F", $2}' ${f_raw}.tmp > "$f_raw"
rm "$f_raw".tmp

cat "${sample}"_all_GAGC.txt "${sample}"_all_GCTC.txt >> "${sample}"_all_GAGCGCTC.txt

cutoff=14
for motif in GAGC GCTC ; do
  awk -v c=$cutoff '$4>=c{print $0}' "$sample"_pileup_dep0."$motif".combine > "$sample"_pileup_dep"$cutoff"."$motif".combine
done

cat "$sample"_pileup_dep"$cutoff".GAGC.combine "$sample"_pileup_dep"$cutoff".GCTC.combine >> "$sample"_pileup_dep"$cutoff".GAGCGCTC.combine

# convert pos to motif.gff
# modified motif
motif=GAGCGCTC
f_pileup="$sample"_pileup_dep"$cutoff"."$motif".combine
fout="$sample"_pileup_dep"$cutoff"_"$motif"_geneClass.gff

sh pileup_to_gffClass.sh $intergene_gff $f_pileup $fout

# all motifs
f_raw="${sample}"_all_"$motif".txt
fout="$sample"_all_"$motif"_geneClass.gff

sh pileup_to_gffClass.sh $intergene_gff $f_raw $fout


# summary gene class.gff
f_motifgff="$sample"_pileup_dep"$cutoff"_"$motif"_geneClass.gff # motif sites in gff format.
sh gffClass_summary_2024.sh $intergene_gff $f_motifgff

f_motifgff="$sample"_all_"$motif"_geneClass.gff # motif sites in gff format.
sh gffClass_summary_2024.sh $intergene_gff $f_motifgff 0

# check numbers.
awk -F '\t' '$6=="tRNA"{print $0}' Lanchnosp_all_GAGCGCTC_geneClass.gff | wc -l
