# Disrupted hybrid gene expression in dwarf hamsters and house mice

This repository is for data processing and analysis related to Hunnicutt et al. 2024, "Different complex regulatory phenotypes underlie hybrid male sterility in divergent rodent crosses" with Colin Callahan, Sara Keeble, Emily C. Moore, Jeffrey M. Good, and Erica L. Larson.

Below, are example commands that were used to take raw RNAseq fastq data to  input invariant + variant VCF files for pixy. All analyses were conducted on the University of Denver Research Data Analysis Cluster (RDAC) running Red Hat Enterprise Linux Server release 7.9 (Maipo). Reference hamster genome is available on NCBI (*Phodopus sungorus*: GCA_023856395.1), and mouse reference was GRCm38.p6.
# Read processing and mapping
## Read trimming and QC (same as for count file generation)
1. FastQC to check quality of raw reads (Illumina NovaSeq 6000 - 150bp PE)
```{bash}
module load tools/fastqc/0.11.9 
mkdir /fastqc_untrimmed/
mkdir /fastqc_untrimmed/sample_fqc-out/

gunzip sample_1.fq.gz 
gunzip sample_2.fq.gz 

fastqc -t 4 -o /fastqc_untrimmed/sample_fqc-out -f fastq sample_1.fq
fastqc -t 4 -o /fastqc_untrimmed/sample_fqc-out -f fastq sample_2.fq

gzip sample_1.fq
gzip sample_2.fq
```

2. Trimmomatic to remove adapters and remove low quality bases/reads
```{bash}
module load tools/trimmomatic/0.39
mkdir trimm-logs
mkdir post-trimm/
mkdir post-trimm/paired/
mkdir post-trimm/unpaired/

#run trimmomatic
java -jar /cm/shared/tools/trimmomatic/0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 \
-trimlog trimm-logs/sample_trimm.log \
sample_1.fq sample_2.fq \
post-trimm/paired/sample_1.paired.fq.gz \
post-trimm/unpaired/sample_1.unpaired.fq.gz \
post-trimm/paired/sample_2.paired.fq.gz  \
post-trimm/unpaired/sample_2.unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36
```

3. FastQC to check quality of trimmed reads
```{bash}
module load tools/fastqc/0.11.9 

mkdir post-trimm/fastqc/sample_fqc-out/
gunzip sample_1.paired.fq.gz
gunzip sample_2.paired.fq.gz

fastqc -t 4 -o post-trimm/fastqc/sample_fqc-out/ \
-f fastq post-trimm/paired/sample_1.paired.fq
fastqc -t 4 -o post-trimm/fastqc/sample_fqc-out/ \
-f fastq post-trimm/paired/sample_2.paired.fq

gzip post-trimm/paired/sample_1.paired.fq
gzip post-trimm/paired/sample_2.paired.fq
```

## Create invariant + variant sites files for hamsters and house mice

4. Map RNASeq reads to reference genomes with hisat2 then sort reads
```{bash}
module load tools/hisat2/2.2.0 
module load tools/samtools/1.10 

#align reads to reference genome and adds read group tag
hisat2 -p 8 -k 500 ref_genome/ref_genome.fna \
--rg-id sample --rg SM:sample -S \
-1 post-trimm/paired/sample_1.paired.fq \
-2 post-trimm/paired/sample_2.paired.fq \
-S post-trimm/post-hisat/example_ref.sam

#convert to bam then sort
samtools view -@ 24 -S -b sample_ref.sam
samtools sort -@ 16 -o sample_ref_sorted.sam -O sample_ref.bam
```

5. Mark Duplicates with Picard
```{bash}
module load tools/picard/2.25.7
cd post-trimm/post-hisat/picard_ref_MD_tmp
mkdir post-trimm/post-hisat/picard_ref_MD_tmp

##marking duplicates with picard
java -jar picard.jar MarkDuplicates \
      TMP_DIR=post-trimm/post-hisat/picard_ref_MD_tmp \
      I=sample_ref.bam\
      O=sample_ref_marked_duplicates.bam \
      M=sample_ref_marked_duplicates_metrics.txt

echo "Optical duplicates marked for sample"
date
```

6. Additional processing needed before HaplotypeCaller for RNASeq data because of splicing (SplitNCigar)
```{bash}
module load tools/gatk/4.2.5.0

gatk SplitNCigarReads \
      -R ref_genome/ref_genome.fna \
      -I sample_ref_marked_duplicates.bam \
      -O sample_ref.split.bam \
      --tmp-dir post-trimm/post-hisat/split_CIGAR_tmpDir
```

7. HaplotypeCaller for each species
```{bash}
module load tools/samtools/1.10 
module load tools/gatk/4.2.5.0

cd post-trimm/post-hisat/
mkdir post-haplotypeCaller

#index bams first
samtools index sample_ref.split.bam

#run HaplotypeCaller
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ref_genome/ref_genome.fna \
   -I sample_ref.split.bam \
   -O post-haplotypeCaller/sample_ref.g.vcf.gz \
   -ERC GVCF
```

8. Make companion files needed for GATK
```{bash}
cd post-trimm/post-hisat/post-haplotypeCaller
#ref
rm ref_RNASeq.sample_map
for q in `ls *_ref.g.vcf.gz`;
do sample=`echo $q | sed 's#_ref.g.vcf.gz##g'`;
file=`echo $q`;
echo $sample '../'$file | tr ' ' '\t' >> ref_RNASeq.sample_map;
done

#need a file with lengths for chroms?
cd ref_genome/
rm ref_genome.length
cat ref_genome.fna.fai | cut -f1-2 | tr '\t' ',' | grep "chr" >> ref_genome.length
```

9. Run GenomicsDBImport for each chromosome for pcam_pseudogenome
```{bash}
module load tools/gatk/4.2.5.0

#make sure these directories exist before running genomicsDBimport
#all GATK steps must by run chromosome by chromosome or you cannot generate the invariant steps file
mkdir post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/
mkdir post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/tmp/


gatk --java-options "-Xmx20g -Xms20g" \
   GenomicsDBImport \
   --genomicsdb-workspace-path ref_chrName/ \
   -L chrName \
   --batch-size 4 \
   --sample-name-map post-trimm/post-hisat/post-haplotypeCaller/ref_RNASeq.sample_map \
   --tmp-dir post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/tmp/ \
   --reader-threads 4
```

10. Run genotypeGVCFs for each chromosome for psun_pseudogenome; allow invariant sites to be called
```{bash}
module load tools/gatk/4.2.5.0
#make sure you create this directory
mkdir post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/post_genotypeGVCFs

cd post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/post_genotypeGVCFs

gatk --java-options "-Xmx20g" GenotypeGVCFs \
   -R ref_genome/ref_genome.fna \
   -V gendb://../ref_chrName/ \
	-L chrName \
   -O ref_RNASeq_ref_chrName_genotypeGVCFs.vcf.gz \
   -all-sites
```

11. Make a reference sequence dictionary prior to filtering
```{bash}
module load tools/picard/2.25.7

cd ref_genome/

java -jar /cm/shared/tools/picard/2.25.7/picard.jar CreateSequenceDictionary \
      R=ref_genome.fna
```

12. Hard filter vcfs with VariantFiltration first (marks variants) then selectVariants (selectsVariants)
```{bash}
module load tools/gatk/4.2.5.0
cd post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/post_genotypeGVCFs
module load tools/gatk/4.2.5.0

#filter variants
gatk VariantFiltration \
   -R ref_genome/ref_genome.fna \
   -V ref_RNASeq_ref_chrName_genotypeGVCFs.vcf.gz \
   -O ref_RNASeq_ref_chrName_filtered.vcf.gz \
   --mask-extension 5 \
   --filter-name "QDFilter" \
   --filter-expression "QD < 2.0" \
   --filter-name "FSFilter" \
   --filter-expression "FS > 60.0" \
   --filter-name "MQFilter" \
   --filter-expression "MQ < 40.0" \
   --filter-name "QUALFilter" \
   --filter-expression "QUAL < 30.0" \
   --filter-name "DP_FilterMin" \
   --filter-expression "DP < 5"

#select variants pcam
 gatk SelectVariants --variant ref_RNASeq_ref_chrName_filtered.vcf.gz \
 --exclude-filtered \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --output RNASeq_ref_chrName_selected.vcf.gz
```

13. Merge chromosomes together after GenotypeGVCFs
```{bash}
module load tools/bcftools/1.10.2 

cd post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/post_genotypeGVCFs

ls *filtered*vcf.gz > ref_RNASeq_ref_allSites_intervals_genotypeGVCFs.list

bcftools concat \
     -f ref_RNASeq_ref_allSites_intervals_genotypeGVCFs.list \
    -o ref_RNASeq_ref_allSites_mergeIntervals_genotypeGVCFs.vcf.gz \
    -O z
```

14. Calculating average depth metrics for use in determining soft filters to eliminate multimapped reads; max coverage allowed was around 2.5x the average
```{bash}
#genome wide coverage for bams

module load tools/samtools/1.17
cd post-trimm/post-hisat/

for q in `ls *split.bam`;
do echo $q;
id=`echo $q | sed 's#_ref.split.bam##g'`;
samtools depth $id'_ref.split.bam' > $id'_ref.split.depth'
echo "generated depth file";
echo $id'_ref.split.depth' >> ref.split.depth
cat $id'_ref.split.depth' | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> ref.split.depth
echo "calculated depth";
done
```

15. Soft-filtering with VCFtools (min and max DP were system dependent; 5 and 60 for hamsters; 5 and 45 for mice); various missingness_thresholds were also tested (Supplemental Figure S3)
```{bash}
cd post-trimm/post-hisat/post-haplotypeCaller/post_genomicsDBimport/post_genotypeGVCFs/
module load tools/vcftools/0.1.17 

vcftools --vcf ref_RNASeq_ref_allSites_mergeIntervals_genotypeGVCFs.vcf.gz \
		--max-missing missingness_threshold \
		--minDP 5 \
		--maxDP 60 \
		--recode \
		--out ref_RNASeq_ref_allSites_finalFilter_maxMissing_missingness_threshold


#gzip then index vcf files
module load tools/htslib/1.13  
module load tools/bcftools/1.10.2cd 
 
bgzip ref_RNASeq_ref_allSites_finalFilter_maxMissing_missingness_threshold.recode.vcf;  
tabix ref_RNASeq_ref_allSites_finalFilter_maxMissing_missingness_threshold.recode.vcf.gz;  
```

## Pixy

16. Create a pixy conda environment; This was run locally on a 2021 MacBook Pro with an Apple M1 Max chip with 64 GB of RAM; Using Conda v. 24.1.2 and python version : 3.11.7.final.0
```{bash}
conda create -n pixy
conda activate pixy

conda install -c conda-forge pixy
conda install -c bioconda htslib

conda deactivate
```

17. Running pixy! This runs using two nested loops that loop through all the invariant files for different missingness thresholds then runs on each chromosome separately. Window size was 10,000bp and 4 cores used
```{bash}
for q in `ls ref_RNASeq_ref_allSites_finalFilter_maxMissing_*gz`;
do echo $q;
output_name=`echo $q | cut -f2 -d"/" | sed 's#.recode.vcf.gz##g'`;
mkdir $output_name;
for chromosome in `cat ref_genome.length | cut -f1 -d","`;
do echo $chromosome;
pixy --stats pi fst dxy \
--vcf $q \
--populations ref_RNASeq_ref_allSites_popfile.txt \
--window_size 10000 \
--n_cores 4 \
--output_folder $output_name \
--output_prefix $chromosome \
--chromosomes $chromosome >> $output_name'/'$chromosome'.log';
done
done
```

18. After Pixy runs, all subsequent analyses took place in R. See accompanying Pixy_calculations.Rmd script.