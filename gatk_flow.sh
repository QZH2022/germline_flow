ID=$1
FASTQ_R1=$2
FASTQ_R2=$3

fastqName=`basename $FASTQ_R1`
prefix=${fastqName%%.*}

home="$ID/res"

REF="Homo_sapiens_assembly38.fasta"
KNOW="/NAS/sl/hg38/"
thread_num=40
dbsnp="dbsnp_146.hg38.vcf.gz"
golden_indel="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KG_phase1_snv="1000G_phase1.snps.high_confidence.hg38.vcf.gz"

mkdir res && mkdir res/coverage && mkdir res/vcf_stats
bts=`date +%s`
echo "stat_start ${prefix} $bts" >> ./time_stat.txt

#######fastp########
fastp="docker run --rm -v $ID:/var/data -v $home:/var/data/output cdgcwgs:3.0 fastp"
time $fastp -i $FASTQ_R1 -I $FASTQ_R2 -o output/${prefix}_clean.R1.fq.gz -O output/${prefix}_clean.R2.fq.gz --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cut_tail --cut_front -l 45 -n 10

#########bwa#######
cd $home
RG="@RG\\tID:group1\\tSM:${prefix}\\tPL:illumina\\tLB:CasCADE\\tPU:unit1"
bwa="docker run --rm  -v $home:/var/data -v $KNOW:/var/data/InputDIR cdgcwgs:3.0 bwa"
time $bwa mem -t 80 -R $RG \
/var/data/InputDIR/$REF \
${prefix}_clean.R1.fq.gz \
${prefix}_clean.R2.fq.gz \
> ${prefix}.sorting
##########samtools#####
samtools="docker run --rm -v $home:/var/data cdgcwgs:3.0 samtools"
time $samtools sort \
--threads 80 \
-O bam ${prefix}.sorting \
> ${prefix}.sorted.bam

#####picard###
picard="docker run --rm -v $home:/var/data -v $KNOW:/var/data/InputDIR
cdgcwgs:3.0 picard"
time $picard AddOrReplaceReadGroups \
I=${prefix}.sorted.bam \
O=${prefix}.sorted.RG.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=${prefix}

picard="docker run --rm -v $home:/var/data -v $KNOW:/var/data/InputDIR
cdgcwgs:3.0 picard"
time $picard \
MarkDuplicates \
INPUT=${prefix}.sorted.RG.bam \
OUTPUT=${prefix}.dedup.bam \
METRICS_${prefix}=${prefix}.dedup.metrics
time $picard \
BuildBamIndex \
INPUT=${prefix}.dedup.bam \
OUTPUT=${prefix}.dedup.bai


############gatk#########
gatk="docker run --rm -v $home:/var/data -v $KNOW/:/var/data/InputDIR
cdgcwgs:3.0 /root/miniconda2/bin/gatk"
time $gatk --java-options "-Xmx60g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
BaseRecalibrator \
-R /var/data/InputDIR/Homo_sapiens_assembly38.fasta \
-I ${prefix}.dedup.bam \
--known-sites /var/data/InputDIR/$dbsnp \
--known-sites /var/data/InputDIR/$golden_indel \
--known-sites /var/data/InputDIR/$KG_phase1_snv \
-O ${prefix}.recal_data.table

time $gatk --java-options "-Xmx60g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
ApplyBQSR \
-I ${prefix}.dedup.bam \
-R  /var/data/InputDIR/Homo_sapiens_assembly38.fasta  \
-bqsr ${prefix}.recal_data.table \
-O ${prefix}.recal.bam


gatk="docker run --rm -v $home:/var/data -v $KNOW/:/var/data/InputDIR
cdgcwgs:3.0 /root/miniconda2/bin/gatk"
time $gatk --java-options "-Xmx60g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
HaplotypeCaller \
-R /var/data/InputDIR/Homo_sapiens_assembly38.fasta \
-I ${prefix}.recal.bam \
-ERC GVCF \
--dbsnp /var/data/InputDIR/$dbsnp \
-O ${prefix}.raw.snps.indels.g.vcf \
--max-reads-per-alignment-start 0

gatk="docker run --rm -v $home:/var/data -v $KNOW/:/var/data/InputDIR
cdgcwgs:3.0 /root/miniconda2/bin/gatk"
time $gatk --java-options "-Xmx60g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
GenotypeGVCFs \
-R /var/data/InputDIR/Homo_sapiens_assembly38.fasta \
-V ${prefix}.raw.snps.indels.g.vcf \
--dbsnp /var/data/InputDIR/$dbsnp \
-O ${prefix}.raw.snps.indels.vcf

ets=`date +%s`
time=`expr $ets - $bts`
echo "stat_1 ${prefix} $time" >> ./time_stat.txt
	
picard="docker run --rm -v $home:/var/data -v $KNOW:/var/data/InputDIR
cdgcwgs:3.0 picard"
time $picard \
CollectWgsMetrics \
I=${prefix}.dedup.bam \
O= coverage/${prefix}.CollectWgsMetrics.txt \
R=/var/data/InputDIR/Homo_sapiens_assembly38.fasta

picard="docker run --rm -v $home:/var/data -v $KNOW:/var/data/InputDIR
cdgcwgs:3.0 picard"
time $picard \
CollectMultipleMetrics \
I=${prefix}.dedup.bam \
O=coverage/${prefix}.CollectMultipleMetrics \
R=/var/data/InputDIR/Homo_sapiens_assembly38.fasta \
PROGRAM=null \
PROGRAM=CollectAlignmentSummaryMetrics \
PROGRAM=CollectInsertSizeMetrics \
PROGRAM=QualityScoreDistribution \
PROGRAM=MeanQualityByCycle \
PROGRAM=CollectBaseDistributionByCycle \
PROGRAM=CollectGcBiasMetrics \
PROGRAM=CollectSequencingArtifactMetrics \
PROGRAM=CollectQualityYieldMetrics \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=temp

bgzip="docker run --rm -v $home:/var/data cdgcwgs:3.0 bgzip"
bcftools="docker run -v $home:/var/data cdgcwgs:3.0 bcftools"
time $bgzip -@ 40 \
-f ${prefix}.raw.snps.indels.vcf
time $bcftools index \
--threads 40 \
-f \
-t ${prefix}.raw.snps.indels.vcf.gz
time $bcftools stats \
${prefix}.raw.snps.indels.vcf.gz \
> vcf_stats/${prefix}.vcf.stats

ets1=`date +%s`
time=`expr $ets1 - $bts`
echo "stat_2 ${prefix} $time" >> ./time_stat.txt