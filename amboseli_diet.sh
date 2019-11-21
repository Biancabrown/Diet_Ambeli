#!/bin/bash

#SBATCH -J met
#SBATCH -o qiime2_o
#SBATCH -e qiime2_e
#SBATCH -n 16
#SBATCH --mem=100GB
#SBATCH -t 24:00:00


#########################################################
#########################################################
#########################################################

####################
#LOAD OSCAR MODULES#
####################

module load anaconda/4.6.14

##
conda activate 

####Get Corn Genome
source activate /users/bbrown3/data/bbrown3/miniconda/envs/ncbi


ncbi-genome-download --genus "Zea mays" plant

module load blast/2.9.0+
export BLASTDB=/users/bbrown3/scratch/
echo $BLASTDB
update_blastdb.pl taxdb
gunzip -cd taxdb.tar.gz | (cd $BLASTDB; tar xvf - )



makeblastdb -in /users/bbrown3/data/bbrown3/amboseli/corn/refseq/chloroplast_maize.fasta -title chloroplast_maize -out /users/bbrown3/data/bbrown3/amboseli/corn/refseq/db_chloroplast_maize -parse_seqids -dbtype nucl


source activate /users/bbrown3/data/bbrown3/miniconda/envs/magicblast
#source activate /users/bbrown3/data/bbrown3/miniconda/envs/metaphlan2

magicblast -query SRR1747065.1_1.fasta -db rbcl.fasta -out rbcl.sam

magicblast -query SRR1747065.1_1.fasta -db rbcl.fasta -outfmt tabular -out rbcl.tab



##Samtools to transfer files sam to bam file for downstream analysis

#!/bin/bash

#SBATCH -J sam
#SBATCH -o sam_o
#SBATCH -e sam_e
#SBATCH -n 16
#SBATCH --mem=100GB
#SBATCH -t 24:00:00

#########################################################
#########################################################
#########################################################

####################
#LOAD OSCAR MODULES#
####################
 module load samtools

 samtools view -S -b corn_output > corn_output.bam

 samtools sort corn_output.bam -o corn_output.sorted.bam

 samtools index corn_output.sorted.bam

source activate /users/bbrown3/data/bbrown3/miniconda/envs/ncbi-genome
ncbi-genome-download --taxid NC_001666 plant




#header
awk 'sub(/^>/, "")' rbcl.fasta > rbcl.txt

sed 's/\s.*$//' rbcl.txt > rbcl_accesion.

grep -c "^>" SRR1747065.1_1.fasta

#reduce samples

#identifying word limit.
#I'm not sure if this best way to do this. 

module load blast/2.9.0+
export BLASTDB=/users/bbrown3/scratch/
echo $BLASTDB
update_blastdb.pl taxdb
gunzip -cd taxdb.tar.gz | (cd $BLASTDB; tar xvf - )

#module load anaconda/4.6.14

##

#conda activate 

#source activate /users/bbrown3/data/bbrown3/miniconda/envs/blast

#blastn -db trnl-f.fas -query SRR1747065.1_1.fasta  -word_size 28 -dust no -perc_identity 90 -num_threads 8 -out blastn_stringent_90.txt -outfmt "6 qseqid sseqid evalue pident length"

blastn -db rbcl.fasta -query SRR1747065.1_1.fasta -word_size 28 -dust no -perc_identity 90  -num_threads 8 -out blastn_rbcl_90_3_taxa.txt -outfmt "6 qseqid sseqid pident length mismatch qstart sstart evalue bitscore staxids sscinames sblastnames"

#Create out put file that can be transferred to R
awk -F'[;\t]' '!seen[$1,$2]++' blastn_rbcl_taxa.txt | awk '{print $1 "\t" $2}' > krona_blastn_rbcl_taxa.txt

