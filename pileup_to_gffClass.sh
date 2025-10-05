#!/bin/bash
gff=$1 # intregene.gff file.
fpos=$2 # pileup txt file.
fout=$3 # output.
gsource=$4 # re-format gid name.
replace=$5

> $fout
# find the line in gff that start <= pos <= end.
awk '{print $1,$2,$4}' $fpos | while read scaffold pos dep; do
  fold=$(echo ${scaffold} | sed "s/^${gsource}/${replace}/")
  LC_ALL=C grep --no-group-separator "^$fold[[:blank:]]" $gff | awk -F '\t' -v p=$pos -v f=$fold -v d=$dep '
  BEGIN { tag=1 ;OFS="\t" }
  { if ($4<=p&&$5>=p) { print f"\t"p"\t"d"\t"$0; tag=2 } }
  END { if ( tag == 1 ) {print "not_found"} }' >> $fout
done
