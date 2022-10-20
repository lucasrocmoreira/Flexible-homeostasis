#!/bin/bash
source /broad/software/scripts/useuse
reuse UGER
reuse Samtools
reuse Java-1.8

# Software paths
fastqc="/seq/vgb/software/FastQC/"
trimmomatic="/seq/vgb/software/Trimmomatic-0.39/"
hisat2="/seq/vgb/software/hisat2/hisat2-2.1.0/"
stringtie="/seq/vgb/software/stringtie/stringtie-1.3.4d.Linux_x86_64/"

# Reference genome
reference="/seq/vgb/lmoreira/reference_genomes/Cer_sim"
annotation="/seq/vgb/lmoreira/reference_genomes/Ceratotherium_simum.gff3"

# If the input is R1, what is the path for R2?
read_path=`dirname $1`
sample_name=`basename $1 | cut -d '_' -f1`
read1=$1
read2=$read_path/$sample_name"_R2_001.fastq.gz"

# File directories 
trimmed_reads="/seq/vgb/flexible_homeostasis/rna_seq/rhino/trimmed_reads"
bamDir="/seq/vgb/flexible_homeostasis/rna_seq/rhino/BAM_files"
stringtieDir="/seq/vgb/flexible_homeostasis/rna_seq/rhino/stringtie"

echo $read1 $read2

###################################
#echo FastQC
###################################

#$fastqc/fastqc $read1
#$fastqc/fastqc $read2

###################################
echo Trimmomatic
###################################

java -jar $trimmomatic/trimmomatic-0.39.jar PE $read1 $read2 \
$trimmed_reads/$sample_name.R1_paired.fq.gz \
$trimmed_reads/$sample_name.R1_unpaired.fq.gz \
$trimmed_reads/$sample_name.R2_paired.fq.gz \
$trimmed_reads/$sample_name.R2_unpaired.fq.gz \
-threads 16 \
ILLUMINACLIP:/seq/vgb/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:true

###################################
echo HiSat2
###################################

$hisat2/hisat2 -p 8 --dta -x $reference -1 $trimmed_reads/$sample_name.R1_paired.fq.gz -2 $trimmed_reads/$sample_name.R2_paired.fq.gz -U $trimmed_reads/$sample_name.R1_unpaired.fq.gz,$trimmed_reads/$sample_name.R2_unpaired.fq.gz | samtools sort -@ 8 -T /seq/vgb/flexible_homeostasis/temp/ -o $bamDir/$sample_name.sorted.bam

###################################
echo Stringtie
###################################

$stringtie/stringtie -o $stringtieDir$id/$sample_name"_stringtie.gtf" -A $stringtieDir$id/$sample_name"_gene_abundances.tsv" -G $annotation -e -v -p 8 $bamDir/$sample_name.sorted.bam

# the option -e measures expression of genes only present in the annotation
# the option -x ignores certain chromosomes (e.g. mt)
