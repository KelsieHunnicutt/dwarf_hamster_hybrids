# Disrupted hybrid gene expression in dwarf hamsters and house mice

This repository is for data processing and analysis related to Hunnicutt et al. 2024, "Different complex regulatory phenotypes underlie hybrid male sterility in divergent rodent crosses" with Colin Callahan, Sara Keeble, Emily C. Moore, Jeffrey M. Good, and Erica L. Larson.

Below, are example commands that were used to take raw RNAseq fastq data to count data. All analyses were conducted on the University of Denver Research Data Analysis Cluster (RDAC) running Red Hat Enterprise Linux Server release 7.9 (Maipo). Reference hamster genome is available on NCBI (*Phodopus sungorus*: GCA_023856395.1).
# Read processing and mapping
## Read trimming and QC
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

## Create *P. sungorus* and a *P. campbelli* pseudogenomes to remove bias associated with mapping hybrids to one reference genome

4. Remove minor scaffolds from NCBI genome then index reference genome with bwa-mem2
```{bash}
mkdir phodopus_genome/

echo CM043852.1 JAJQIY010000284.1 JAJQIY010007291.1 CM043853.1 JAJQIY010004270.1 JAJQIY010000046.1 JAJQIY010006224.1 JAJQIY010007057.1 CM043854.1 JAJQIY010003985.1 JAJQIY010007356.1 CM043855.1 JAJQIY010000470.1 JAJQIY010002681.1 JAJQIY010003188.1 JAJQIY010003390.1 JAJQIY010003611.1 JAJQIY010005987.1 JAJQIY010007692.1 CM043856.1 JAJQIY010004751.1 JAJQIY010000586.1 JAJQIY010007049.1 CM043857.1 CM043858.1 JAJQIY010006557.1 CM043859.1 JAJQIY010000219.1 JAJQIY010005426.1 CM043860.1 JAJQIY010002690.1 CM043861.1 JAJQIY010004677.1 CM043862.1 CM043863.1 JAJQIY010006728.1 CM043864.1 JAJQIY010002969.1 > scaffolds_to_keep.txt

#retrieve chromosomes
for scaffold in `cat scaffolds_to_keep.txt`;
do echo $scaffold
grep -A1 scaffold phodopus_genome/GCA_023856395.1_Psun_UM_1.1_genomic.fna >> phodopus_genome/psun_genome.fna;
done

module load tools/bwa-mem2/2.0 
cd phodopus_genome/

bwa-mem2 index psun_genome.fna
samtools faidx psun_genome.fna
```

5. Download Pcam WGS and Psun RNASeq data from NCBI
```{bash}
module load tools/sra-tools/2.10.8

mkdir WGS-data/
cd WGS-data/

#retrieve files from ncbi then convert .sra to .fastq
prefetch SRR17223279 --max-size 100GB
prefetch SRR17223284 --max-size 100GB

fasterq-dump SRR17223279.sra --split-files
fasterq-dump SRR17223284.sra --split-files
```

6. Map Pcam WGS (NCBI SRA: SRR17223279) reads to Psun genome with bwa-mem2 then sort reads
```{bash}
module load tools/bwa-mem2/2.0 
module load tools/samtools/1.10

mkdir WGS-data/
cd WGS-data/

#align reads to reference genome
bwa-mem2 mem -t 16 \
-K 100000000 -p -v 3 -t 16 -Y \
phodopus-genomes/psun_genome.fna \
SRR17223279_1.fastq SRR17223279_2.fastq \
-o pcam_aligned_to_psun_NCBI.bam

#sort bam
samtools sort -@ 16 -o pcam_aligned_to_psun_NCBI_sorted.bam -O bam pcam_aligned_to_psun_NCBI.bam
```

7. Map Psun RNASeq (NCBI SRA: SRR17223284/SSSS_231-6M_WT) reads to Psun genome with hisat2 then sort reads
```{bash}
module load tools/hisat2/2.2.0 
module load tools/samtools/1.10 

cd WGS-data/

#align reads to reference genome and adds read group tag
hisat2 -p 8 -k 500 phodopus_genome/psun_genome.fna \
--rg-id SSSS_231-6M_WT --rg SM:SSSS_231-6M -S \
-1 post-trimm/paired/SSSS_231-6M_WT_1.paired.fq \
-2 post-trimm/paired/SSSS_231-6M_WT_2.paired.fq \
-S post-trimm/post-hisat/ref/SSSS_231-6M_WT_ref.sam

#convert to bam then sort
samtools view -@ 24 -S -b SSSS_231-6M_WT_ref.sam
samtools sort -@ 16 -o SSSS_231-6M_WT_ref_sorted.sam -O SSSS_231-6M_WT_ref.sam
```

8. Mark Duplicates with Picard
```{bash}
module load tools/picard/2.25.7
cd WGS-data/

#pcam
java -jar /cm/shared/tools/picard/2.25.7/picard.jar MarkDuplicates \
      TMP_DIR=pcam_dup_tmp \
      I=pcam_aligned_to_psun_NCBI_sorted.bam \
      O=pcam_aligned_to_psun_marked_duplicates.bam \
      M=pcam_aligned_to_psun_marked_dup_metrics.txt

#psun
java -jar /cm/shared/tools/picard/2.25.7/picard.jar MarkDuplicates \
      TMP_DIR=mark_duplicates_SSSS_WT_tmpDir \
      I=SSSS_231-6M_WT_ref_sorted.sam \
      O=psun_aligned_to_psun_marked_duplicates.bam \
      M=psun_aligned_to_psun_marked_dup_metrics.txt
```

9. Additional processing needed before HaplotypeCaller for RNASeq data because of splicing (SplitNCigar)
```{bash}
module load tools/gatk/4.2.5.0

gatk SplitNCigarReads \

      -R phodopus_genome/psun_genome.fna \
      -I psun_aligned_to_psun_marked_duplicates.bam \
      -O psun_aligned_to_psun.split.bam \
      --tmp-dir WGS-data/split_CIGAR_SSSS_WT_tmpDir
```

10. HaplotypeCaller for each species
```{bash}
module load tools/samtools/1.10 
module load tools/gatk/4.2.5.0

cd WGS-data/

#index bams first
samtools index psun_aligned_to_psun.split.bam
samtools index pcam_aligned_to_psun_marked_duplicates.bam

#run haplotype caller
#pcam
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R phodopus_genome/psun_genome.fna \
   -I pcam_aligned_to_psun.marked_duplicates.sorted.bam \
   -O pcam_aligned_to_psun.g.vcf.gz \
   -ERC GVCF
   
#psun
gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R phodopus_genome/psun_genome.fna \
   -I psun_aligned_to_psun.split.bam \
   -O psun_aligned_to_psun.g.vcf.gz \
   -ERC GVCF
```

11. Run genotypeGVCFs for each chromosome for pcam_pseudogenome
```{bash}
module load tools/gatk/4.2.5.0

mkdir WGS-data/genotypeGVCFs_byChr/
cd WGS-data/genotypeGVCFs_byChr/

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R phodopus_genome/psun_genome.fna \
   -L numberedChr \
   -V ../pcam_aligned_to_psun.g.vcf.gz \
   -O pcam_aligned_to_psun_numberedChr_genotypedGVCF.vcf.gz

#make numberedChr specific
for q in `cat phodopus_genome/psun_genome.fna.fai | cut -f1`;
do echo $q;
sed 's#numberedChr#'$q'#g' name_of_file_containing_the_GenotypeGVCF_pcam_script.txt > gatk_genotypeGVCFs_pcam_$q'.submit';
done
```

12. Run genotypeGVCFs for each chromosome for psun_pseudogenome
```{bash}
module load tools/gatk/4.2.5.0

mkdir WGS-data/genotypeGVCFs_byChr/
cd WGS-data/genotypeGVCFs_byChr/

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R phodopus_genome/psun_genome.fna \
   -L numberedChr \
   -V ../psun_aligned_to_psun.g.vcf.gz \
   -O psun_aligned_to_psun_numberedChr_genotypedGVCF.vcf.gz

#make numberedChr specific
for q in `cat phodopus_genome/psun_genome.fna.fai | cut -f1`;
do echo $q;
sed 's#numberedChr#'$q'#g' name_of_file_containing_the_GenotypeGVCF_script.txt > gatk_genotypeGVCFs_psun_$q'.submit';
done
```

13. Make a reference sequence dictionary prior to filtering
```{bash}
module load tools/picard/2.25.7

cd phodopus_genome/

java -jar /cm/shared/tools/picard/2.25.7/picard.jar CreateSequenceDictionary \
      R=psun_genome.fna
```

14. Hard filter vcfs with VariantFiltration first (marks variants) then selectVariants (selectsVariants)
```{bash}
module load tools/gatk/4.2.5.0
cd WGS-data/genotypeGVCFs_byChr/

#filter variants pcam
gatk VariantFiltration \
   -R phodopus_genome/psun_genome.fna \
   -V pcam_aligned_to_psun_numberedChr_genotypedGVCF.vcf.gz \
   -O pcam_aligned_to_psun_numberedChr_filtered.vcf.gz \
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
   --filter-expression "DP < 10" \
   --filter-name "DP_FilterMax" \
   --filter-expression "DP > 150"

#filter variants psun
gatk VariantFiltration \
   -R phodopus_genome/psun_genome.fna \
   -V psun_aligned_to_psun_numberedChr_genotypedGVCF.vcf.gz \
   -O psun_aligned_to_psun_numberedChr_filtered.vcf.gz \
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
   --filter-expression "DP < 10" \
   --filter-name "DP_FilterMax" \
   --filter-expression "DP > 150"

#select variants pcam
 gatk SelectVariants --variant pcam_aligned_to_psun_numberedChr_filtered.vcf.gz \
 --exclude-filtered \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --output pcam_aligned_to_psun_numberedChr_selected.vcf.gz

#select variants psun
 gatk SelectVariants --variant psun_aligned_to_psun_numberedChr_filtered.vcf.gz \
 --exclude-filtered \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --output psun_aligned_to_psun_numberedChr_selected.vcf.gz
```

15. Incorporate SNPs back into reference genome to make pseudogenomes for each chromosome
```{bash}
module load tools/gatk/4.2.5.0
cd WGS-data/genotypeGVCFs_byChr/

#pcam
gatk FastaAlternateReferenceMaker \
 -R phodopus_genome/psun_genome.fna \
 -O phodopus_genome/pcam_numberedChr_pseudogenome.fasta \
 -L numberedChr \
 -V pcam_aligned_to_psun_numberedChr_selected.vcf.gz 

#pcam
gatk FastaAlternateReferenceMaker \
 -R phodopus_genome/psun_genome.fna \
 -O phodopus_genome/psun_numberedChr_pseudogenome.fasta \
 -L numberedChr \
 -V psun_aligned_to_psun_numberedChr_selected.vcf.gz 
```

16. Concatenate  chromosome pseudogenomes
```{bash}
#pcam
cd phodopus_genome/
for q in `ls pcam_*numberedChr*_pseudogenome.fasta`; 
do echo $q; 
cat $q >> pcam_pseudogenome.fna;
done

#psun
cd phodopus_genome/
for q in `ls pcam_*numberedChr*_pseudogenome.fasta`; 
do echo $q; 
cat $q >> pcam_pseudogenome.fna;
done

#fix pseudogenome headers so that both pseudogenomes match
sed -i -e 's#>1 #>#' pcam_to_NCBI_allChr_pseudogenome.fasta
sed -i -e 's#:.*##' pcam_to_NCBI_allChr_pseudogenome.fasta
```

17. Generate a mod file using lapels and modtools from pseudogenome for use in modtools pipeline later
```{bash}
module load tools/picard/2.25.7
module load apps/python/2.7.18
source /cm/shared/apps/python/2.7.18/packages/lapels/bin/activate
cd phodopus_genome/

#first need to create metadata files from reference genome for both pseudogenomes
get_refmeta -o Pcam.1.meta Pcam.1 psun_genome.fasta
get_refmeta -o Psun.1.meta Psun.1 psun_genome.fasta

#concatenate chromosome level vcfs
cd WGS-data/genotypeGVCFs_byChr/
java -jar picard.jar MergeVcfs \
          I=list_of_selected_pcam_vcfs_by_chr.list \
          O=pcam_allChr_merged_selected.vcf

java -jar picard.jar MergeVcfs \
          I=list_of_selected_psun_vcfs_by_chr.list \
          O=psun_allChr_merged_selected.vcf
          
#vcf2mod
cd WGS-data/

#pcam
vcf2mod -v -o pcam_pseudogenome.mod phodopus_genome/psun_genome.fna Pcam.1.meta PCam.1 \
WGS-data/genotypeGVCFs_byChr/pcam_allChr_merged_selected.vcf

#psun
vcf2mod -v -o psun_pseudogenome.mod phodopus_genome/psun_genome.fna Psun.1.meta Psun.1 \
WGS-data/genotypeGVCFs_byChr/psun_allChr_merged_selected.vcf
```

## Mapping reads
18. Index the reference genome and the pseudogenomes and create an easier file name alias
```{bash}
module load tools/hisat2/2.2.0 
cd phodopus_genome

hisat2-build -p 24 psun_genome.fna psun_genome.fna
hisat2-build -p 24 psun_genome.fna psun_pseudogenome.fna
hisat2-build -p 24 psun_genome.fna psun_pseudogenome.fna
```

19. Align reads to both pseudogenomes with HiSat2. 
```{bash}
module load tools/hisat2/2.2.0 
module load tools/samtools/1.10 
mkdir post-trimm/post-hisat/psun_ref
mkdir post-trimm/post-hisat/psun_pseudo
mkdir post-trimm/post-hisat/pcam_pseudo

#align to P. sungorus pseudogenome
hisat2 -p 4 -k 250 -x phodopus_genome/psun_genome.fna \
-1 post-trimm/paired/sample_1.paired.fq \
-2 post-trimm/paired/paired/sample_2.paired.fq \
-S post-trimm/post-hisat/psun_ref/sample_psun_ref.sam

#align to P. sungorus pseudogenome
hisat2 -p 4 -k 250 -x phodopus_genome/psun_pseudogenome.fna \
-1 post-trimm/paired/sample_1.paired.fq \
-2 post-trimm/paired/paired/sample_2.paired.fq \
-S post-trimm/post-hisat/psun_pseudo/sample_psun_pseudo.sam

#align to P. campbelli pseudogenome
hisat2 -p 4 -k 250 -x phodopus_genome/pcam_pseudogenome.fna \
-1 post-trimm/paired/sample_1.paired.fq \
-2 post-trimm/paired/paired/sample_2.paired.fq \
-S post-trimm/post-hisat/pcam_pseudo/sample_pcam_pseudo.sam

#convert sams to bams
cd post-trimm/post-hisat/psun_pseudo/
samtools view -S -b sample_psun_ref.sam > sample_psun_ref.bam
rm sample_psun_ref.sam

cd post-trimm/post-hisat/psun_pseudo/
samtools view -S -b sample_psun_pseudo.sam > sample_psun_pseudo.bam
rm sample_psun_pseudo.sam

cd post-trimm/post-hisat/pcam_pseudo/
samtools view -S -b sample_pcam_pseudo.sam > sample_pcam_pseudo.bam
rm sample_pcam_pseudo.sam
```

20. Sort bam files
```{bash}
module load tools/samtools/1.10 

samtools sort -o sample_ref_sorted.bam -O bam -@ 8 sample_psun_ref.bam

samtools sort -o sample_psun_sorted.bam -O bam -@ 8 sample_pcam_pseudo.bam

samtools sort -o sample_pcam_sorted.bam -O bam -@ 8 sample_pcam_pseudo.bam
```

21. Append HI (hit index; distinguishes mappings of multi-mapped reads) to bam files with the hisat2tophat.py script from https://github.com/goodest-goodlab/pseudo-it/blob/master/helper-scripts/hisat2Tophat.py
```{bash}
module load apps/python/3.6

#append HI
python3 hisat2tophat.py sample_ref_sorted.bam sample_ref_postHI.bam
python3 hisat2tophat.py sample_ref_sorted.bam sample_psun_postHI.bam
python3 hisat2tophat.py sample_ref_sorted.bam sample_pcam_postHI.bam
```

22. Run lapels on pcam and psun mapping file to convert back to reference genome coordinates
```{bash}
source /cm/shared/apps/python/2.7.18/packages/lapels/bin/activate

pylapels -p 24 -n -o sample_pcam_post-mod.bam phodopus_genome/pcam_pseudogenome.mod \
sample_pcam_postHI.bam

pylapels -p 24 -n -o sample_psun_post-mod.bam phodopus_genome/psun_pseudogenome.mod \
sample_psun_postHI.bam

deactivate
```

23. Run suspenders to merge pcam and psun alignments and keep best alignments for each read
```{bash}
module load tools/suspenders/0.2.6
module load tools/samtools/1.10 

pysuspenders ./sample_post-sus.bam \
sample_pcam_post-mod.bam  \
sample_psun_post-mod.bam 

#sort
samtools sort -o sample_post-sus_sorted.bam \
-O bam -@ 2 sample_post-sus.bam
```

23. FeatureCounts needs a gtf file instead of gff so convert with gffread
```{bash}
cd phodopus_genome/

gffread psun_annotations.gff -T -o psun_annotations_postGffread.gtf
```

24. Generate gene count data with featureCounts. 
```{bash}
module load tools/subread/2.0.1 

#allow multi-mapping (not main file used in manuscript)
featureCounts -B -C -M -T 4 \
-a phodopus_genome/psun_annotations_postGffread.gtf \
-o Hamsters_psun_14Mar23_multiMapping_featureCounts.csv \
BBBB_153-5M_DIP_post-sus_sorted.bam BBBB_153-5M_LZ_post-sus_sorted.bam BBBB_153-5M_SP_post-sus_sorted.bam BBBB_153-6M_DIP_post-sus_sorted.bam BBBB_153-6M_SP_post-sus_sorted.bam BBBB_154-4M_LZ_post-sus_sorted.bam BBBB_195-4M_DIP_post-sus_sorted.bam BBBB_195-4M_LZ_post-sus_sorted.bam BBBB_195-4M_RS_post-sus_sorted.bam BBBB_195-4M_SP_post-sus_sorted.bam BBBB_197-3M_DIP_post-sus_sorted.bam BBBB_197-3M_RS_post-sus_sorted.bam BBBB_197-4M_DIP_post-sus_sorted.bam BBBB_197-4M_LZ_post-sus_sorted.bam BBBB_197-4M_RS_post-sus_sorted.bam BBBB_198-7M_LZ_post-sus_sorted.bam BBBB_198-7M_RS_post-sus_sorted.bam BBBB_198-7M_SP_post-sus_sorted.bam BBBB_280-4M_WT_post-sus_sorted.bam BBBB_281C-3M_WT_post-sus_sorted.bam BBBB_283B-4M_WT_post-sus_sorted.bam BBSS_100-2M_DIP_post-sus_sorted.bam BBSS_100-2M_LZ_post-sus_sorted.bam BBSS_100-2M_SP_post-sus_sorted.bam BBSS_7-1M_WT_post-sus_sorted.bam BBSS_90-2M_DIP_post-sus_sorted.bam BBSS_90-2M_LZ_post-sus_sorted.bam BBSS_90-2M_SP_post-sus_sorted.bam BBSS_90-3M_SP_post-sus_sorted.bam BBSS_95-3M_DIP_post-sus_sorted.bam BBSS_95-3M_LZ_post-sus_sorted.bam BBSS_99-3M_DIP_post-sus_sorted.bam BBSS_99-3M_LZ_post-sus_sorted.bam BBSS_99-3M_SP_post-sus_sorted.bam SSSS_138-5M_DIP_post-sus_sorted.bam SSSS_138-5M_LZ_post-sus_sorted.bam SSSS_155-5M_DIP_post-sus_sorted.bam SSSS_155-5M_LZ_post-sus_sorted.bam SSSS_155-5M_SP_post-sus_sorted.bam SSSS_155-6M_DIP_post-sus_sorted.bam SSSS_155-6M_LZ_post-sus_sorted.bam SSSS_155-6M_SP_post-sus_sorted.bam SSSS_182-5M_DIP_post-sus_sorted.bam SSSS_182-5M_LZ_post-sus_sorted.bam SSSS_182-5M_RS_post-sus_sorted.bam SSSS_188-4M_DIP_post-sus_sorted.bam SSSS_188-4M_RS_post-sus_sorted.bam SSSS_188-4M_SP_post-sus_sorted.bam SSSS_188-5M_RS_post-sus_sorted.bam SSSS_188-5M_SP_post-sus_sorted.bam SSSS_231-6M_WT_post-sus_sorted.bam SSSS_233A-7M_WT_post-sus_sorted.bam SSSS_233B-4M_WT_post-sus_sorted.bam```

#do not allow multi-mapping (main file used in manuscript)
featureCounts -B -C -T 12 -a /data/larsonlab/publicResources/genomes/hamster/ncbi/dwarf_hamster_annotation_with_NCBI_names_postGffread.gtf \

-o Hamsters_psun_14Mar23_singleMapping_featureCounts.csv \

BBBB_153-5M_DIP_post-sus_sorted.bam BBBB_153-5M_LZ_post-sus_sorted.bam BBBB_153-5M_SP_post-sus_sorted.bam BBBB_153-6M_DIP_post-sus_sorted.bam BBBB_153-6M_SP_post-sus_sorted.bam BBBB_154-4M_LZ_post-sus_sorted.bam BBBB_195-4M_DIP_post-sus_sorted.bam BBBB_195-4M_LZ_post-sus_sorted.bam BBBB_195-4M_RS_post-sus_sorted.bam BBBB_195-4M_SP_post-sus_sorted.bam BBBB_197-3M_DIP_post-sus_sorted.bam BBBB_197-3M_RS_post-sus_sorted.bam BBBB_197-4M_DIP_post-sus_sorted.bam BBBB_197-4M_LZ_post-sus_sorted.bam BBBB_197-4M_RS_post-sus_sorted.bam BBBB_198-7M_LZ_post-sus_sorted.bam BBBB_198-7M_RS_post-sus_sorted.bam BBBB_198-7M_SP_post-sus_sorted.bam BBBB_280-4M_WT_post-sus_sorted.bam BBBB_281C-3M_WT_post-sus_sorted.bam BBBB_283B-4M_WT_post-sus_sorted.bam BBSS_100-2M_DIP_post-sus_sorted.bam BBSS_100-2M_LZ_post-sus_sorted.bam BBSS_100-2M_SP_post-sus_sorted.bam BBSS_7-1M_WT_post-sus_sorted.bam BBSS_90-2M_DIP_post-sus_sorted.bam BBSS_90-2M_LZ_post-sus_sorted.bam BBSS_90-2M_SP_post-sus_sorted.bam BBSS_90-3M_SP_post-sus_sorted.bam BBSS_95-3M_DIP_post-sus_sorted.bam BBSS_95-3M_LZ_post-sus_sorted.bam BBSS_99-3M_DIP_post-sus_sorted.bam BBSS_99-3M_LZ_post-sus_sorted.bam BBSS_99-3M_SP_post-sus_sorted.bam SSSS_138-5M_DIP_post-sus_sorted.bam SSSS_138-5M_LZ_post-sus_sorted.bam SSSS_155-5M_DIP_post-sus_sorted.bam SSSS_155-5M_LZ_post-sus_sorted.bam SSSS_155-5M_SP_post-sus_sorted.bam SSSS_155-6M_DIP_post-sus_sorted.bam SSSS_155-6M_LZ_post-sus_sorted.bam SSSS_155-6M_SP_post-sus_sorted.bam SSSS_182-5M_DIP_post-sus_sorted.bam SSSS_182-5M_LZ_post-sus_sorted.bam SSSS_182-5M_RS_post-sus_sorted.bam SSSS_188-4M_DIP_post-sus_sorted.bam SSSS_188-4M_RS_post-sus_sorted.bam SSSS_188-4M_SP_post-sus_sorted.bam SSSS_188-5M_RS_post-sus_sorted.bam SSSS_188-5M_SP_post-sus_sorted.bam SSSS_231-6M_WT_post-sus_sorted.bam SSSS_233A-7M_WT_post-sus_sorted.bam SSSS_233B-4M_WT_post-sus_sorted.bam
```

25. After countfile generation, all subsequent analyses took place in R. See accompanying .Rmd scripts for differential expression analyses.