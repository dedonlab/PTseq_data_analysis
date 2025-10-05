#!/bin/bash
lst_g=$1
N=$(wc -l $lst_g | cut -d' ' -f1)
dir_p=$(cd "$(dirname "$2")"; pwd -P)/$(basename "$2")
dir_meme=${dir_p}/memetmp
mkdir $dir_meme || true
sample=$3
#sample=$(basename "$3" | awk -F '[_.]' '{print $1}')

flank=6
flnk2=$(($flank+$flank+1))

#### meme required modules
perlbrew init
perlbrew use perl-5.18.2

module add c3ddb/automake/1.15
module add c3ddb/autoconf/2.69
module add c3ddb/python/3.5.2
module add c3ddb/zlib/1.2.11
module add c3ddb/jdk/1.8.0-101
module add c3ddb/zlib/1.2.11 ; module add c3ddb/xz/5.2.3 ; module add c3ddb/lzma/4.32.7
#export TZ='EST' date
export PATH=/scratch/users/yfyuan/bin/ghostscript-9.52/bin:$PATH
export PATH=/scratch/users/yfyuan/bin/meme-5.3.3/bin:/scratch/users/yfyuan/bin/meme-5.3.3/libexec/meme-5.3.3:$PATH

#### MAIN ####
#<<'####'
for file in $(cat $lst_g) ; do
  #  # remove short squeuences, meme require >= 8
  #  # dreme require same lengtn
    /scratch/users/yfyuan/bin/seqkit/seqkit seq -m 8 $file > ${file}min8
done

####

#<<'####'
for file in $(cat $lst_g) ; do
  # meme.
  # -cefrac 0.8 not used with classic mode
  meme -dna -objfun classic -nmotifs 10 -mod zoops -evt 0.05 -time 3000 -minw 3 -maxw 5 -markov_order 0 -nostatus -oc ${dir_meme} ${file}min8
  mv "$dir_meme"/meme.txt "$dir_p"/meme_"$(basename ${file}|cut -d'.' -f1)".txt
  >&2 echo "$file"
done
