#!/bin/bash
set -euo pipefail

# === Logging setup ===
log_file="pipeline_$(date +%F_%T).log"
exec > >(tee -i "$log_file")
exec 2>&1

echo "=== Starting SNPs Processing Pipeline ==="
echo "Timestamp: $(date)"

# ===Config===
input_vcf="./TotalRawSNPs.vcf"
REF="../reference.fasta"
PICARD_JAR="/tools/picard.jar"
HTSLIB="/tools/htslib-1.21/bin"
GATK="/tools/gatk-4.6.1.0/gatk"
BCFtool="/tools/bcftools-1.21/bin/bcftools"
VCFtool="/tools/vcftools/bin/vcftools"
PLINK="/tools/plink_1.9/plink"
ref_dict="./reference.dict"
export PATH=/tools/jdk-23.0.1/bin:$PATH
#export PATH=/tools/jre1.8.0_311/bin:$PATH

# ===Tool Availability====
echo " Checking required tools..."

for tool in java "$HTSLIB/bgzip" "$HTSLIB/tabix" "$GATK" "$VCFtool" "$BCFtool" "$PLINK"; do
    if command -v "$tool" &>/dev/null || [[ -x "$tool" ]]; then
        echo "Found executable: $tool"
    else
        echo "Error: Required tool not found or not executable -> $tool"
        exit 1
    fi
done

# ===check vcf file in dir===
if [[ ! -f "$input_vcf" ]]; then
	echo "Error: file '$input_vcf' not found"
	exit 1
fi

# ===prepare working dir===
mkdir -p ./vcf_file
cd vcf_file

# ===move input into vcf_file===
echo "[step-1] moving '$input_vcf' into vcf_files"
ln -s ../TotalRawSNPs.vcf .

# ===Compress and index vcf===
echo "[Step-2] Compressing and indexing input VCF..."

$HTSLIB/bgzip -c -@ 16 ./TotalRawSNPs.vcf > TotalRawSNPs.vcf.gz

$HTSLIB/tabix TotalRawSNPs.vcf.gz
#echo "compressing and indexing done"

# ===CreateSequenceDictionary (Picard)=== 
echo "[Step-3] Creating refrence dictionary by picard ..."

java -jar $PICARD_JAR CreateSequenceDictionary R=$REF O=$ref_dict

# ===UpdateVcfSequenceDictionary (Picard) ### contig in vcf headers===
echo "[Step-4] UpdateVcfSequenceDictionary by picard..."

java -jar $PICARD_JAR UpdateVcfSequenceDictionary INPUT=TotalRawSNPs.vcf.gz OUTPUT=TotalRawSNPs_r1.vcf.gz SEQUENCE_DICTIONARY=$ref_dict

# ===index file needed in gatk===
echo "[Step-5] Indexing TotalRawSNPs_r1.vcf.gz for gatk..." 

/$HTSLIB/tabix TotalRawSNPs_r1.vcf.gz

# ===extract Raw SNPs===
echo "[Step-6] Extracting SNPs from VCF using GATK..."

$GATK SelectVariants -R $REF -V TotalRawSNPs_r1.vcf.gz --select-type-to-include SNP -O RawSNPsOnly.vcf.gz --exclude-non-variants --remove-unused-alternates >> RawSNPsOnly.logs 2>&1 

# ===extract BIALLELIC (Step - 1)===
echo "[Step-7]  Extracting biallelic SNPs..."

$GATK SelectVariants -R $REF -V RawSNPsOnly.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O SNPsOnly_r1.vcf.gz >> SNPsOnly_r1.log 2>&1

# ===Filter SNPs based on allele frequency and count===
echo "[Step-8] Filtering SNPs using min and max allele frequency at 0.0001 and 0.9999 ..."

$VCFtool --gzvcf SNPsOnly_r1.vcf.gz --non-ref-af 0.0001 --max-non-ref-af 0.9999 --mac 1 --recode --recode-INFO-all --stdout | $HTSLIB/bgzip -c -@6 > SNPsOnly_r2.vcf.gz

# ===remove loci with Qual < 20===
echo "[Step-9] Filter SNPs based on quality > 20 ..."

$BCFtool filter -e 'QUAL<20' -o SNPsOnly_r3.vcf.gz -O z --threads 8 SNPsOnly_r2.vcf.gz 

# ===indexing file needed for gatk===
echo "[Step-10] Indexing of SNPsOnly_r3.vcf.gz..."

$HTSLIB/tabix SNPsOnly_r3.vcf.gz

# ===Filter for DP<3===
echo "[Step-11] Filtering genotypes based on their depth ..."

$GATK VariantFiltration -R $REF -V SNPsOnly_r3.vcf.gz -O SNPsOnly_r4.vcf.gz -G-filter "DP<3" --genotype-filter-name "dp_lt3" --set-filtered-genotype-to-no-call >> SNPsOnly_r4.log 2>&1

echo "[Step-12] Excluding Filter, non variant i.e. monomorphic and unused alternates  ..."

$GATK SelectVariants -R $REF -V SNPsOnly_r4.vcf.gz -O SNPsOnly_r5.vcf.gz --exclude-filtered --exclude-non-variants --remove-unused-alternates >> SNPsOnly_r5.log 2>&1

# ===remove missing===
echo "[Step-13] Allow 20% missing data, if more remove the loci  ..."

$VCFtool --gzvcf SNPsOnly_r5.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | $HTSLIB/bgzip -c -@ 6 > SNPsOnly_r6.vcf.gz

# ===Minor allele frequency===
echo "[Step-14] Calculating the minor allele frequency (MAF)..."

$PLINK --allow-extra-chr --double-id --freq --vcf SNPsOnly_r6.vcf.gz --out SNPsOnly_r6_maf_stat

# ===filter snps for MAF===
echo "[Step-15-a] Keep only SNPs with MAF 0.01 ..."
$VCFtool --gzvcf SNPsOnly_r6.vcf.gz --maf 0.01 --recode --recode-INFO-all --stdout | $HTSLIB/bgzip -c -@6 > SNPsOnly_r7_01.vcf.gz

echo "[Step-15-b] Keep only SNPs with MAF 0.001 ..."
$VCFtool --gzvcf SNPsOnly_r6.vcf.gz --maf 0.001 --recode --recode-INFO-all --stdout | $HTSLIB/bgzip -c -@6 > SNPsOnly_r7_001.vcf.gz

echo "[Step-15-c] Keep only SNPs with MAF 0.05 ..."
$VCFtool --gzvcf SNPsOnly_r6.vcf.gz --maf 0.05 --recode --recode-INFO-all --stdout | $HTSLIB/bgzip -c -@6 > SNPsOnly_r7_05.vcf.gz

echo "[Step-16] Calculating the stats of all vcf.gz files "
for file in ./*.vcf.gz; do
	name=$(basename "$file" .vcf.gz)
	$BCFtool stats -s - "$file" > "${name}.stat"
done

#=========================================Choose best maf vcf file for filter the heterozygous================================================
