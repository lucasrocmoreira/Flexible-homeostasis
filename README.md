## RNA-seq analysis

Pipeline for the analysis of RNA-seq data for the *Flexible Homeostasis Project*

***

1. Initiate by checking the quality of the sequenced reads using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
fastqc $read1
fastqc $read2
```

2) Trim reads to remove Illumina adaptors using [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```
java -jar trimmomatic-0.39.jar PE $read1 $read2 \
$sample_name.R1_paired.fq.gz \
$sample_name.R1_unpaired.fq.gz \
$sample_name.R2_paired.fq.gz \
$sample_name.R2_unpaired.fq.gz \
-threads 16 \
ILLUMINACLIP:/seq/vgb/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:8:true
```
> * ILLUMINACLIP: cut adapter and other illumina-specific sequences from the read.
> * fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters.
> * seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed [2]
> * palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment [30]
> * simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read. [10]
> * minimum adapter length in palindrome mode [2]
> * keepBothReads [true]

3) Create reference index for [HiSat2](http://daehwankimlab.github.io/hisat2/) alignment

```
hisat2-build $ref.fa $name -p 16
```

> * NOTE: make sure the fasta file for genome reference is not .gz.

We ran on the Broad cluster using:

```
qsub -pe smp 16 -binding linear:16 -l h_vmem=32g -l h_rt=12:00:00 -b y -N hisat2-build -cwd -j y -V /seq/vgb/software/hisat2/hisat2-2.1.0/hisat2-build /home/unix/lmoreira/vgb/reference_genomes/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa /seq/vgb/lmoreira/reference_genomes/mRatBN7.2 -p 16
```

4) Use [HiSat2](http://daehwankimlab.github.io/hisat2/) to align reads to reference genome. 

```
$hisat2/hisat2 -p 8 --dta -x $reference \
-1 $sample_name.R1_paired.fq.gz \
-2 $sample_name.R2_paired.fq.gz \
-U $sample_name.R1_unpaired.fq.gz,$sample_name.R2_unpaired.fq.gz | samtools sort -@ 8 -T /seq/vgb/flexible_homeostasis/temp/ -o $bamDir/$sample_name.sorted.bam
```

5) Quantify transcript abundance using [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml).

```
$stringtie -o $sample_name"_stringtie.gtf" -A $sample_name"_gene_abundances.tsv" -G $annotation.gtf -e -v -p 8 $sample_name.sorted.bam
```
> * NOTE: make sure the gtf file for genome annotation is not .gz.
> * -e measures expression of genes only present in the annotation
> * -v turns on verbose mode

***

* We wrote a simple [script](https://github.com/lucasrocmoreira/Flexible-homeostasis/blob/main/RNA_seq_alignment.sh) to run this pipeline at the Broad server for several samples.

```
for i in vgb/flexible_homeostasis/rna_seq/rat/RN*_R1_001.fastq.gz; 
do name=`basename $i | cut -d '_' -f1`; ./run_job.sh $name vgb/flexible_homeostasis/rna_seq/scripts/RNA_seq_alignment.sh $i; 
done
```
or (when in interactive mode [ish])
```
ish -N RNA-seq -pe smp 16 -binding linear:16 -l h_vmem=32g -l h_rt=12:00:00

for i in vgb/flexible_homeostasis/rna_seq/rat/RN*_R1_001.fastq.gz; 
do vgb/flexible_homeostasis/rna_seq/scripts/RNA_seq_alignment.sh $i; 
done
```

***

### Differentially expession analysis

Based on [this tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).
