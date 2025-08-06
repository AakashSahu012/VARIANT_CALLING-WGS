# Variant Calling Pipeline (Bash-based)

This repository contains a fully automated variant calling pipeline written in Bash. It processes raw whole-genome resequencing (WGS) paired-end FASTQ files to generate high-confidence SNP and InDel calls. The pipeline is designed for scalability, reproducibility, and ease of use on local or HPC systems.

---

## 📁 Directory Structure

project-root/
├── raw_data/ # Input FASTQ files (R1/R2)
├── ref/ # Reference genome (FASTA)
├── scripts/ # Bash scripts for each pipeline step
├── results/ # Final outputs (BAM, VCF, stats)
├── .gitignore
├── README.md
└── run_pipeline.sh # Main executable script


---

## ⚙️ Tools and Dependencies

Ensure the following tools are installed and added to your `PATH`:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)
- [`BWA`](http://bio-bwa.sourceforge.net/)
- [`SAMtools`](http://www.htslib.org/)
- [`BCFtools`](http://samtools.github.io/bcftools/)
- [`GATK`](https://gatk.broadinstitute.org/)
- [`VCFtools`](https://vcftools.github.io/index.html)
- `bash`, `awk`, `sed`, `gzip`

---

## 🔁 Pipeline Workflow

1. **Quality Control**
   - Run `FastQC` on raw FASTQ files.

2. **Trimming**
   - Trim adapters and low-quality bases using `Trimmomatic`.

3. **Alignment**
   - Index the reference genome using `BWA`.
   - Align reads using `bwa mem`.

4. **Post-Processing**
   - Convert SAM to BAM, sort, and index with `SAMtools`.

5. **Variant Calling**
   - Use `BCFtools` or `GATK HaplotypeCaller` for SNP/Indel calling.

6. **Filtering**
   - Apply quality filters using `GATK`, `VCFtools`, or `BCFtools`.

7. **Downstream**
   - Generate filtered VCF and summary statistics.
   - Prepare data for tools like `TASSEL` for trait association.

---

## 🚀 Running the Pipeline

```bash
chmod +x ngs_pipeline.sh
./ngs_pipeline.sh

