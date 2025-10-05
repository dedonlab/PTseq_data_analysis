#!/bin/bash

gene_class_summary(){
  gene_class=$1
  intergeneic_gff=$2 # intergene annoted gff.
  fout=$3 # pileup gff.
  cutoff=$4

  # fout columns: 1.PT_scaffold, 2.PT_pos, 3.scaffold, 4.platform, 5.feature, 6.start, 7.end, 8.., 9.strand, 10.0, 11.description.

  # total number of genes in the class
  ngenegenome=$(awk -v g=$gene_class -F '\t' '$3==g {print $0}' $intergeneic_gff | wc -l )

  # split modification/motif in each gene class
  awk -F '\t' -v c=$cutoff -v g=$gene_class 'BEGIN { OFS = "\t"} $3>=c&&$6==g {print $0}' $fout > ${fout}.${gene_class}

  if [ -s ${fout}.${gene_class} ] ; then
    # number of unique genes modified.
    ngenemod=$(awk '{print $4,$5,$6,$7,$8,$9,$10,$11,$12}' ${fout}.${gene_class} | sort -u | wc -l)

    # proportion that modified
    propmod=$(echo "scale=4; ${ngenemod}/${ngenegenome}" | bc)

    # total PT gene length.
    totalgenelength=$(awk -F '\t' '{print $4,$5,$6,$7,$8,$9,$10,$11,$12}' ${fout}.${gene_class} | sort -u | awk '{print $5-$4}'  | paste -sd+ | bc)

    # average length per gene
    meanlength=$(echo "scale=4; ${totalgenelength}/${ngenemod}" | bc)

    # number of pt sites in the class
    nsites=$(wc -l ${fout}.${gene_class} | cut -d' ' -f1)

    # average number of pt sites per gene
    meansites=$(echo "scale=4; ${nsites}/${ngenemod}" | bc)

    # average distance of sites in the genes
    meandistance=$(python /home/yfyuan/data/scripts/motifscreen/average_dist_in_gene.py ${fout}.${gene_class})
  else
    ngenemod='0'
    propmod='0'
    totalgenelength='0'
    meanlength='0'
    nsites='0'
    meansites='0'
    meandistance='0'
  fi

  echo -e "${gene_class}\t${ngenemod}\t${propmod}\t${totalgenelength}\t${meanlength}\t${nsites}\t${meansites}\t${meandistance}"
}

f_intergenic_gff=$1 # intergene annoted gff.
f_pileup_gff=$2 # pileup gff.
pileup_cutoff=$3

echo -e "Gene_Class\t#gene_modified\tfraction_of_mod_gene\ttotal_gene_length\tmean_length\t#sites\tmean_sites\tmean_distance"
GENECLASS=('intergene' 'tRNA' 'rRNA' 'CDS' 'tmRNA' 'CRISPR' 'repeat_region' 'assembly_gap')
for geneclass in ${GENECLASS[*]} ; do
  gene_class_summary  $geneclass $f_intergenic_gff $f_pileup_gff $pileup_cutoff
done
