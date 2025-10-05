#!/bin/bash

# dependence
# seqkit
# meme v5.3.3
# meme required modules
# perl recommend version 5.18.2
# automake recommend version 1.15
# autoconf recommend version 2.69
# python recommend version 3.5.2
# zlib recommend version 1.2.11
# jdk recommend version 1.8.0
# zlib recommend version 1.2.11
# xz recommend version 5.2.3
# lzma recommend version 4.32.7
# export TZ='EST' date
# ghostscript recommend version 9.52

mkdir meme || true

fas=$1
fas8=${file%.fasta}_min8.fasta

flank=6
flnk2=$(($flank+$flank+1))

# remove short sequences <=8.
# remove short squeuences, meme require >= 8
# dreme require same length

seqkit seq -m 8 ${fas} > ${fas8}

# meme.
# -cefrac 0.8 not used with classic mode
meme -dna -objfun classic -nmotifs 10 -mod zoops -evt 0.05 -time 3000 -minw 3 -maxw 5 -markov_order 0 -nostatus -oc meme ${fas8}

# rename file
mv meme/meme.txt meme/"$(basename ${fas}|cut -d'.' -f1)".txt
