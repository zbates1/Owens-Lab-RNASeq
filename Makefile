# Directories
RAW_DIR := raw_data/
QC_DIR := qc_results/
TRIMMED_DIR := trimmed/
ALIGNED_DIR := aligned/
COUNTS_DIR := counts/

# Files
REF_GENOME := reference_genome.fasta
ANNOTATION := annotation.gtf
INDEX := reference_index

# Rules
all: counts

qc: $(RAW_DIR)/*.fastq
	mkdir -p $(QC_DIR)
	fastqc -o $(QC_DIR) $^

trim: $(RAW_DIR)/*.fastq
	mkdir -p $(TRIMMED_DIR)
	trim_galore --paired -o $(TRIMMED_DIR) $(RAW_DIR)/sample_R1.fastq $(RAW_DIR)/sample_R2.fastq

index:
	hisat2-build $(REF_GENOME) $(INDEX)

align: trim
	mkdir -p $(ALIGNED_DIR)
	hisat2 -x $(INDEX) -1 $(TRIMMED_DIR)/sample_R1_val_1.fq -2 $(TRIMMED_DIR)/sample_R2_val_2.fq -S $(ALIGNED_DIR)/sample.sam
	samtools view -bS $(ALIGNED_DIR)/sample.sam | samtools sort -o $(ALIGNED_DIR)/sample_sorted.bam
	samtools index $(ALIGNED_DIR)/sample_sorted.bam

counts: align
	mkdir -p $(COUNTS_DIR)
	featureCounts -T 4 -a $(ANNOTATION) -o $(COUNTS_DIR)/counts.txt $(ALIGNED_DIR)/*.bam
