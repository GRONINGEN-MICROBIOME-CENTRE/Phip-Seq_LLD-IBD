#!/bin/sh
#SBATCH  --job-name=blast.job
#SBATCH --output=blast.out
#SBATCH --error=blast.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1

ml Anaconda3
#conda create --prefix ./Blast_Env 

source activate  ./Blast_Env 
#conda install -c bioconda blast

INPUT=TGG.fa
INPUT=Probes_arno.fa
INPUT=VFDB_annotated_sb.fa

OUT=Blast_results.txt
OUT=Blast_results_arno.txt
OUT=Blast_results_VFDB.txt

#makeblastdb -in Probes.fa -parse_seqids -dbtype prot
blastp -query $INPUT  -db Probes.fa -out $OUT -outfmt 7 #-evalue 1e-5




