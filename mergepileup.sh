sample=Lanchnosp

# columns: 1scaffold, 2pos, 3cov, 4depR2, 5dep/cov, 6seq.
txtF=${sample}_pileup_dep0_F.pos.txt
txtR=${sample}_pileup_dep0_R.pos.txt

#cout G*AGC +/- 2 nt.
# xxxGAGCTCxxxxx GCTC -> GCTC
# xxxxGAGCTCxxxx GAGC,GCTC -> GCTC
# xxxxxGAGCTCxxx GAGC,GCTC -> GAGC
# xxxxxxGAGCTCxx GAGC,GCTC -> GAGC
# xxxxxxxGAGCTCx GAGC -> GAGC

motif=GAGC
#F
grep " ....${motif}" $txtF | grep -v " ....${motif}TC" > "$txtF".GAGC.tmp
grep " .....${motif}" $txtF >> "$txtF".GAGC.tmp
grep " ......${motif}" $txtF >> "$txtF".GAGC.tmp
grep " .......${motif}" $txtF >> "$txtF".GAGC.tmp
grep " ........${motif}" $txtF >> "$txtF".GAGC.tmp

# remove duplicates.
awk '!a[$0]++' "$txtF".GAGC.tmp > "$txtF".GAGC.tmp.dedup

#R
grep " ....${motif}" $txtR | grep -v " ....${motif}TC" > "$txtR".GAGC.tmp
grep " .....${motif}" $txtR >> "$txtR".GAGC.tmp
grep " ......${motif}" $txtR >> "$txtR".GAGC.tmp
grep " .......${motif}" $txtR >> "$txtR".GAGC.tmp
grep " ........${motif}" $txtR >> "$txtR".GAGC.tmp

# remove duplicates.
awk '!a[$0]++' "$txtR".GAGC.tmp > "$txtR".GAGC.tmp.dedup

#cout G*CTC.
motif=GCTC
#F
grep " ....${motif}" $txtF > "$txtF".GCTC.tmp
grep " .....${motif}" $txtF | grep -v " .....${motif}TC" >> "$txtF".GCTC.tmp
grep " ......${motif}" $txtF | grep -v " ......${motif}TC" >> "$txtF".GCTC.tmp
grep " .......${motif}" $txtF >> "$txtF".GCTC.tmp
grep " ........${motif}" $txtF >> "$txtF".GCTC.tmp

# remove duplicates.
awk '!a[$0]++' "$txtF".GCTC.tmp > "$txtF".GCTC.tmp.dedup

#R
grep " ....${motif}" $txtR > "$txtR".GCTC.tmp
grep " .....${motif}" $txtR | grep -v " .....${motif}TC"  >> "$txtR".GCTC.tmp
grep " ......${motif}" $txtR | grep -v " ......${motif}TC" >> "$txtR".GCTC.tmp
grep " .......${motif}" $txtR >> "$txtR".GCTC.tmp
grep " ........${motif}" $txtR >> "$txtR".GCTC.tmp

# remove duplicates.
awk '!a[$0]++' "$txtR".GCTC.tmp > "$txtR".GCTC.tmp.dedup

# remove duplicates between GAGC and GCTC.
awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtF".GCTC.tmp.dedup "$txtF".GAGC.tmp.dedup > "$txtF".GAGC.tmp.dedup2
awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtR".GCTC.tmp.dedup "$txtR".GAGC.tmp.dedup > "$txtR".GAGC.tmp.dedup2

mv "$txtF".GAGC.tmp.dedup2 "$txtF".GAGC.tmp.dedup
mv "$txtR".GAGC.tmp.dedup2 "$txtR".GAGC.tmp.dedup


# get rest seq
awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtF".GAGC.tmp.dedup "$txtF" > "$txtF".other.tmp
awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtF".GCTC.tmp.dedup "$txtF".other.tmp > "$txtF".other.tmp.tmp
mv "$txtF".other.tmp.tmp "$txtF".other.tmp

awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtR".GAGC.tmp.dedup "$txtR" > "$txtR".other.tmp
awk 'FNR==NR {a[$0];next} !($0 in a)' "$txtR".GCTC.tmp.dedup "$txtR".other.tmp > "$txtR".other.tmp.tmp
mv "$txtR".other.tmp.tmp "$txtR".other.tmp

# merge/collapse closed pileups.
#F
for motif in GAGC GCTC ; do
  f_F="$txtF"."$motif".tmp.dedup
  fout="$txtF"."$motif".tmp.dedup.collapse
  # combine N-2,N-1,N,N+1,N+2.
  # F strand 6->5 +1.
  # columns: 1scaffold, 2pos, 3cov, 4depR2, 5dep/cov, 6seq.
  awk -v r_2=^....${motif} -v r_1=^.....${motif} -v r=^......${motif} -v r1=^.......${motif} -v r2=^........${motif} \
  '{if ($6~r) print $0 ; \
  else if ($6~r_2) print $1,($2+2),$3,$4,$5,"NN"$6;\
  else if ($6~r_1) print $1,($2+1),$3,$4,$5,"N"$6;\
  else if ($6~r1) print $1,($2-1),$3,$4,$5,substr($6,2,length($6)-1)"N";\
  else if ($6~r2) print $1,($2-2),$3,$4,$5,substr($6,2,length($6)-1)"NN"}' "$f_F" | sort -k1,2 > "$fout".tmp

  # combine cov and depth at the same pos ($2).
  awk 'FNR==1{a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}\
  FNR>1{if ($2==b) {c=c+$3;d=d+$4;e=e+$5 } else {print a,b,c,d,e,f ; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}}\
  END {if ($2==b) {c=c+$3;d=d+$4;e=e+$5; print a,b,c,d,e,f} else {print a,b,c,d,e,f"\n"$0} }' "$fout".tmp > $fout
  rm ${fout}.tmp
done


#R
for motif in GAGC GCTC ; do
  f_R="$txtR"."$motif".tmp.dedup
  fout="$txtR"."$motif".tmp.dedup.collapse
  # R strand 6->5 "-" 1
  awk -v r_2=^....${motif} -v r_1=^.....${motif} -v r=^......${motif} -v r1=^.......${motif} -v r2=^........${motif} \
  '{if ($6~r) print $0 ; \
  else if ($6~r_2) print $1,($2-2),$3,$4,$5,"NN"$6;\
  else if ($6~r_1) print $1,($2-1),$3,$4,$5,"N"$6;\
  else if ($6~r1) print $1,($2+1),$3,$4,$5,substr($6,2,length($6)-1)"N";\
  else if ($6~r2) print $1,($2+2),$3,$4,$5,substr($6,2,length($6)-1)"NN"}' "$f_R"| sort -k1,2 > "$fout".tmp

  # combine cov and depth at the same pos ($2).
  awk 'FNR==1{a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}\
  FNR>1{if ($2==b) {c=c+$3;d=d+$4;e=e+$5 } else {print a,b,c,d,e,f ; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}}\
  END {if ($2==b) {c=c+$3;d=d+$4;e=e+$5; print a,b,c,d,e,f} else {print a,b,c,d,e,f"\n"$0} }' "$fout".tmp > $fout
  rm ${fout}.tmp
done



# 3. concat F R files
for motif in GAGC GCTC ; do
  f_Fmotif="$txtF"."$motif".tmp.dedup.collapse
  f_Rmotif="$txtR"."$motif".tmp.dedup.collapse
  fout=${sample}_pileup_dep0.${motif}.combine
  # concat F, R txt files.
  awk -v m=$motif 'FNR==NR{print $1,$2,$3,$4,($4/$3),$6," F "m} FNR<NR{print $1,$2,$3,$4,($4/$3),$6," R "m}' $f_Fmotif $f_Rmotif | sort -t ' ' -Vk1,2 > ${fout}
done

# 3. concat other sequces to motif F R files
# combine cov and depth at the pos +/-3 nt.
awk 'FNR==1{a=$1;b=$2;c=$3;d=$4;e=$5;f=$6} \
FNR>1{if ($1==a&&$2<=b+2&&$4>=d) {b=$2;c=c+$3;d=d+$4;e=e+$5;f=$6} \
else if ($1==a&&$2<=b+2&&$4<d) {c=c+$3;d=d+$4;e=e+$5} \
else {print a,b,c,d,e,f ; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}} \
END {if ($1==a&&$2<=b+2&&$4>=d) {b=$2;c=c+$3;d=d+$4;e=e+$5;f=$6;print a,b,c,d,e,f} \
else if ($1==a&&$2<=b+2&&$4<d) {c=c+$3;d=d+$4;e=e+$5; print a,b,c,d,e,f} else {print a,b,c,d,e,f"\n"$0}}' "$txtF".other.tmp > "$txtF".other.collapse


awk 'FNR==1{a=$1;b=$2;c=$3;d=$4;e=$5;f=$6} \
FNR>1{if ($1==a&&$2<=b+2&&$4>=d) {b=$2;c=c+$3;d=d+$4;e=e+$5;f=$6} \
else if ($1==a&&$2<=b+2&&$4<d) {c=c+$3;d=d+$4;e=e+$5} \
else {print a,b,c,d,e,f ; a=$1;b=$2;c=$3;d=$4;e=$5;f=$6}} \
END {if ($1==a&&$2<=b+2&&$4>=d) {b=$2;c=c+$3;d=d+$4;e=e+$5;f=$6;print a,b,c,d,e,f} \
else if ($1==a&&$2<=b+2&&$4<d) {c=c+$3;d=d+$4;e=e+$5; print a,b,c,d,e,f} else {print a,b,c,d,e,f"\n"$0}}' "$txtR".other.tmp > "$txtR".other.collapse


# concat other sequces to motif F R files
awk 'FNR==NR{print $0" F other"} FNR<NR{print $0" R "m}' \
"$txtF".other.tmp \
"$txtR".other.tmp | sort -t ' ' -Vk1,2 > \
${sample}_pileup_dep0.other.combine


rm *tmp
