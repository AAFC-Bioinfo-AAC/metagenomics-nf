# Metagenomic_nf Documentation Index

---

## Description

The Metagenomic_nf pipeline is a Nextflow-based workflow designed for automated metagenomics analyses, including quality filtering, taxonomic classification, functional profiling, assembly, binning, and annotation of metagenome-assembled genomes (MAGs). This documentation provides a structured guide to setting up, running, and customizing the pipeline. Users will find references for dependencies, environment setup, and detailed workflow descriptions to facilitate execution across different computing environments.

---

## Table of Contents

- [Description](#description)
- [Technical Resources](#technical-resources)
  - [Configuration Files](#configuration-files)
  - [Tool Documentation](#tool-documentation)
  - [Environment Setup](#environment-setup)
- [Glossary](#glossary)
  - [Metagenomics and Bioinformatics](#metagenomics-and-bioinformatics)
  - [Software and Tools](#software-and-tools)
- [Additional Documentation Files](#additional-documentation-files)
  - [README.md](#readmemd)
  - [Metagenomic_nf Guide](#metagenomic_nf-guide)

---

## Technical Resources

### Configuration Files
The pipeline relies on configuration files to define parameters, input data locations, and execution settings. Key files include:

- `nextflow.config`: Defines profiles, resource allocations, environment paths, and parameters for pipeline execution.
- `metadata/map_file.tsv`: Provides sample-to-file mapping for raw reads.
- `metadata/coassembly_groups.tsv`: Defines co-assembly sample groupings.
- `metadata/config.json`: Specifies configurations for functional annotation tools like DRAM.
- `parameters.config`: Contains pipeline parameter overrides for customized execution.

### Tool Documentation
Below are key tools and bioinformatics software used within the pipeline:

- [Nextflow](https://www.nextflow.io/)
- [Fastp](https://github.com/OpenGene/fastp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)
- [Kaiju](https://github.com/bioinformatics-centre/kaiju)
- [HUMAnN](https://huttenhower.sph.harvard.edu/humann)
- [MEGAHIT](https://github.com/voutcn/megahit)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- [CheckM2](https://github.com/chklovski/CheckM2)
- [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/)
- [PhyloPhlAn](https://github.com/biobakery/phylophlan)
- [CoverM](https://github.com/wwood/CoverM)
- [DRAM](https://github.com/WrightonLabCSU/DRAM)

### Environment Setup
For detailed environment setup, refer to the [Usage](/README.md#1-dependencies) section in the README. The pipeline is designed to be executed within Conda environments specified in `nextflow.config`.

Prebuilt Conda environments are available for specific computing environments (e.g., AAFC/AAC Biocluster, PHAC/NML Waffles Cluster, and GPSC). The required packages are detailed in [conda_env.md](/docs/conda_env.md).

---

## Glossary

### Metagenomics and Bioinformatics
- **Metagenomics**: The study of genetic material recovered directly from environmental samples.
- **Taxonomic Profiling**: The classification of microbial species within a sample using tools like Kraken2 and Kaiju.
- **Functional Profiling**: The annotation of metabolic pathways and gene functions using tools like HUMAnN and DRAM.
- **Assembly & Binning**: The process of reconstructing genomes from metagenomic sequencing data using assemblers like MEGAHIT and binning tools like MetaBAT2.
- **MAGs (Metagenome-Assembled Genomes)**: Genomic reconstructions obtained from metagenomic data through binning and refinement processes.

### Software and Tools
- **Bowtie2**: A tool for aligning sequencing reads to a reference genome.
- **Kraken2 & Kaiju**: Tools for taxonomic classification of metagenomic reads.
- **MEGAHIT & MetaBAT2**: Assemblers and binning tools for reconstructing genomes from metagenomic data.
- **CheckM2 & GTDB-Tk**: Tools for quality assessment and phylogenetic classification of MAGs.
- **PhyloPhlAn & CoverM**: Tools for phylogenetic profiling and metagenomic read coverage analysis.
- **DRAM**: A tool for annotating microbial metabolism and gene functions.

---

## Additional Documentation Files

### README.md
Refer to [README.md](/README.md) for an overview of the repository, pipeline execution instructions, and environment setup details.

### ABCC_RCBA_Guide
Refer to the [ABCC_RCBA_Guide](https://github.com/AAFC-Bioinformatics/ABCC_RCBA_Guide) for additional context and supplementary materials that align with this and other AAFC projects.

---

For further inquiries or contributions, refer to the [CONTRIBUTING.md](/CONTRIBUTING.md) file.