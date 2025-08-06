#!/usr/bin/bash

#### Create only once ####

# Define directories
input_dir="/home/iiita/ngs_data/data/input"
output_dir="/home/iiita/ngs_data/data/output_dir"
trimmed_dir="$output_dir/trimmed"
aligned_dir="$output_dir/aligned"
results_dir="$output_dir/results"

# Create necessary directories
mkdir -p $output_dir $trimmed_dir $aligned_dir $results_dir

# Define reference genome
reference_genome="/home/iiita/aakash_ref/genome.fasta"

#### Indexing the reference genome ####
## indexing of .fa file for alignment 
bwa index reference_genome

########## END ###########



####### Start Analysis #######

# Step 1: FastQC - Quality check of raw data
echo "Running FastQC on raw data..."
fastqc $input_dir/*.fastq.gz -o $output_dir


##Step 2: Trimmomatic - Trimming low quality reads and adapters
echo "Running Trimmomatic..."
for file in $input_dir/*_R1_001.fastq.gz; do
    base=$(basename "$file" _R1_001.fastq.gz)
    java -jar /usr/bin/trimmomatic-0.39.jar PE \
        "$input_dir/${base}_R1_001.fastq.gz" "$input_dir/${base}_R2_001.fastq.gz" \
        "$trimmed_dir/${base}_R1_paired.fastq.gz" "$trimmed_dir/${base}_R1_unpaired.fastq.gz" \
        "$trimmed_dir/${base}_R2_paired.fastq.gz" "$trimmed_dir/${base}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
done


# Step 3: BWA-MEM - Aligning reads to the reference genome
echo "Running BWA-MEM..."
for file in $trimmed_dir/*_R1_paired.fastq.gz; do
    base=$(basename "$file" _R1_paired.fastq.gz)
    bwa mem "$reference_genome" "$trimmed_dir/${base}_R1_paired.fastq.gz" "$trimmed_dir/${base}_R2_paired.fastq.gz" > "$aligned_dir/${base}.sam"
done


# Step 4: SAMtools - Convert SAM to BAM, sort, and index
echo "Processing SAM files with SAMtools..."
for file in $aligned_dir/*.sam; do
    base=$(basename "$file" .sam)
    # Convert SAM to BAM, sort, and save as .sorted.bam
    samtools view -bS "$file" | samtools sort -o "$aligned_dir/${base}.sorted.bam"
    # Index the sorted BAM file with CSI index
    samtools index -c "$aligned_dir/${base}.sorted.bam"
done


# Step 5: BCFtools - Generate VCF
echo "Running BCFtools mpileup and call..."
for file in $aligned_dir/*.sorted.bam; do
    base=$(basename "$file" .sorted.bam)
   bcftools mpileup -Ou -f "$reference_genome" "$file" | bcftools call -mv -Oz -o "$results_dir/${base}.vcf.gz"
done


# Step 6: Index the VCF files
echo "Indexing VCF files..."
for file in $results_dir/*.vcf.gz; do
    bcftools index "$file"
done

echo "Pipeline completed successfully. Results are in $results_dir"

