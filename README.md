## Copy right and Contact
PI: pcdedon@mit.edu

Author and contact: yfyuan@mit.edu, msdemott@mit.edu

## Edit log
10/04/2025 revised: add test dataset, version number

## Overview
This foler contains bash and python scripts and bash commands used in the manuscripts entitled "Phosphorothioate DNA modification by BREX type 4 systems in the human gut microbiome" and "PT-seq: A method for metagenomic analysis of phosphorothioate epigenetics in complex microbial communities"

It aims to align reads of an example PT-seq dataset, identify read pileups, and extract sequences including 6 flanking nt at the pileup sites. The package utilizes a simply pipeline from trimming to mapping, sorting and sequencing analysis.

## Hardware requirements
The PTseq pileup analysis pipelien requires only a standard computer with enough RAM to support the in-memory operations.

## OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

Linux: Centos 7
macOS: Mojave (10.14.1)

## Installation and dependences
bbmap https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/

fastqc v0.11.8 https://github.com/s-andrews/FastQC

bowtie2 v2.4.5 https://github.com/BenLangmead/bowtie2/releases

samtools v1.19.2 https://github.com/samtools/samtools/releases/

bedtools v2.30.0 https://github.com/arq5x/bedtools2/releases


MEME-suit v5.3.3 #Please follow the instructions of the MEME suite via its website at http://meme-suite.org

Dependences for MEME-suit v5.3.3

perl v5.18.2 https://www.cpan.org/src/

ghostscript v9.52 https://ghostscript.com/releases/

automake v1.15 https://ftp.gnu.org/gnu/automake/

autoconf v2.69 https://www.gnu.org/software/autoconf/

python v3.5.2 https://www.python.org/downloads/release/python-352/

zlib v1.2.11 https://github.com/jrmwng/zlib-1.2.11

jdk v1.8.0-101 https://www.oracle.com/java/technologies/javase/javase8-archive-downloads.html

zlib v1.2.11 https://github.com/jrmwng/zlib-1.2.11

xz v5.2.3 https://github.com/tukaani-project/xz/releases

lzma v4.32.7 https://sourceforge.net/projects/lzma/

## Usage
1. install the dependences
2. Download the scripts and the demo dataset. Place them in the work directory
3. Keep the demo dataset in work/demo/
4. To deplete human sequence contamination (optional), please download the hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz file at https://zenodo.org/records/1208052 and place it in the work/demo folder with the demo reads.
5. Modify trim.sh with ${path_to_your_bbmap}
6. Run the scripts in the following order:

    \# RAM >= 50G is required. For real PT-seq dataset, we recommond thread >= 10.
   
    `bash sh trim.sh demo/demo_1.fastq demo/demo_2.fastq job_demo`

The output files for the next step are trimmed reads: job_demo_R1_final.fastq job_demo_R2_final.fastq. The output files also include intermediate .fq files and QC report files.
   
   `sh main.sh`
   
   `sh pos2seq_R.sh`
   
   `sh mergepileup.sh`
   
   `sh meme.sh`
   
   `sh motif_stat.sh`
   
   `sh pileup_to_gffClass.sh`
   
   `sh summary_geneClass.sh`
   
   `sh gene_class_summary.sh`
