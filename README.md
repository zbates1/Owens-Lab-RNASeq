# Welcome to Owen's Lab RNASeq Pipeline!

## Introduction
This pipeline supports both **bulk-RNASeq** and **single-cell RNASeq (scRNASeq)** analyses using cutting-edge tools. Below, you'll find an overview of the technologies and workflows.

## Tools
Here are the tools used in the pipeline:

| Analysis Type       | Tool             | Documentation Link                                                                                              |
|---------------------|------------------|----------------------------------------------------------------------------------------------------------------|
| Bulk-RNASeq         | PyDESeq2        | [PyDESeq2 Docs](https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_minimal_pydeseq2_pipeline.html)    |
| Single-Cell RNASeq  | SCANPY (TBD)    | [SCANPY Docs]()                                                                                                 |

The single-cell RNASeq analysis leverages the **SCANPY** package, while **PyDESeq2** is a new (Sept 2023) Python alternative to the widely used DESeq2 R library. PyDESeq2 achieves performance on par with DESeq2, making it an excellent option for Python users.

### Comparison of PyDESeq2 and DESeq2
![Comparison of PyDESeq2 and DESeq2](/images/pydeseq2.png)

[Source](https://academic.oup.com/bioinformatics/article/39/9/btad547/7260507)

---

## Workflow
The workflow for **bulk-RNASeq** involves the following steps:

1. Quality Control
2. Alignment
3. Metadata Curation
4. Differential Expression Analysis with PyDESeq2

![Bulk-RNASeq Workflow Diagram](/images/workflow.png)

---

## Data Requirements
The PyDESeq2 library requires:
- **Input counts**: Raw gene expression counts
- **Metadata**: Sample information for experimental design

Prepare these inputs carefully using quality control and alignment tools.

---

## CLI Packages to create input counts and metadata
```bash
sudo apt update
```

```bash
sudo apt install fastqc -y
```

```bash
sudo apt install python3-pip -y
```

```bash
pip3 install cutadapt
```

```bash
sudo apt install trim-galore -y
```

```bash
sudo apt install hisat2 -y
```

```bash
sudo apt install subread -y
```

```bash
sudo apt install samtools -y
```

---

## Downloading Reference Genome and Annoation File
Mus musculus (Mouse, GRCm39):

```bash
wget ftp://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz
```

```bash
# Decompress files
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.109.gtf.gz
```

Prepare these inputs carefully using quality control and alignment tools.

---


## Next Steps
- [ ] Complete documentation links for scRNASeq tools.
- [ ] Add more images or diagrams for the scRNASeq workflow.
