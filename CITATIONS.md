# Citations for Metagenomic_nf Pipeline

This document contains citations for tools, databases, and methodologies used in the **Metagenomic_nf** pipeline. Please ensure proper attribution when using this pipeline in scientific publications.

---

### **Assembly & Binning**
- **MEGAHIT**: Li, D., Liu, C-M., Luo, R., et al. (2015). "MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph." *Bioinformatics*, 31(10), 1674-1676. [DOI:10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)
- **MetaBAT2**: Kang, D. D., Li, F., Kirton, E., et al. (2019). "MetaBAT 2: An adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies." *PeerJ*, 7, e7359. [DOI:10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359)
- **CheckM2**: Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). "CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes." *Genome Research*, 25(7), 1043-1055. [DOI:10.1101/gr.186072.114](https://doi.org/10.1101/gr.186072.114)

### **Taxonomic Classification**
- **Kraken2**: Wood, D. E., Lu, J., & Langmead, B. (2019). "Improved metagenomic analysis with Kraken 2." *Genome Biology*, 20, 257. [DOI:10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)
- **Bracken**: Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). "Bracken: Estimating species abundance in metagenomics data." *PeerJ Computer Science*, 3, e104. [DOI:10.7717/peerj-cs.104](https://doi.org/10.7717/peerj-cs.104)
- **Kaiju**: Menzel, P., Ng, K. L., & Krogh, A. (2016). "Fast and sensitive taxonomic classification for metagenomics with Kaiju." *Nature Communications*, 7, 11257. [DOI:10.1038/ncomms11257](https://doi.org/10.1038/ncomms11257)

### **Functional Annotation**
- **HUMAnN3**: Franzosa, E. A., McIver, L. J., Rahnavard, G., et al. (2018). "HUMAnN 2.0: A computational framework for improving functional profiling of microbial communities." *Nature Methods*, 15(11), 881-888. [DOI:10.1038/s41592-018-0176-y](https://doi.org/10.1038/s41592-018-0176-y)
- **DRAM**: Shaffer, M., Borton, M. A., McGivern, B. B., et al. (2020). "DRAM for distilling microbial metabolism to automate the curation of microbiome function." *Nucleic Acids Research*, 48(16), 8883-8900. [DOI:10.1093/nar/gkaa621](https://doi.org/10.1093/nar/gkaa621)

### **Phylogenetics & Dereplication**
- **GTDB-Tk**: Chaumeil, P. A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). "GTDB-Tk: A toolkit to classify genomes with the Genome Taxonomy Database." *Bioinformatics*, 36(6), 1925-1927. [DOI:10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)
- **PhyloPhlAn**: Segata, N., Börnigen, D., Morgan, X. C., & Huttenhower, C. (2013). "PhyloPhlAn is a new method for improved phylogenetic and taxonomic placement of microbes." *Nature Communications*, 4, 2304. [DOI:10.1038/ncomms3304](https://doi.org/10.1038/ncomms3304)
- **dRep**: Olm, M. R., Brown, C. T., Brooks, B., & Banfield, J. F. (2017). "dRep: A tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes." *The ISME Journal*, 11, 2864-2868. [DOI:10.1038/ismej.2017.126](https://doi.org/10.1038/ismej.2017.126)

### **Quality Control & Read Processing**
- **Fastp**: Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). "fastp: An ultra-fast all-in-one FASTQ preprocessor." *Bioinformatics*, 34(17), i884-i890. [DOI:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- **Bowtie2**: Langmead, B., & Salzberg, S. L. (2012). "Fast gapped-read alignment with Bowtie 2." *Nature Methods*, 9, 357-359. [DOI:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)
- **Bedtools**: Quinlan, A. R., & Hall, I. M. (2010). "BEDTools: A flexible suite of utilities for comparing genomic features." *Bioinformatics*, 26(6), 841-842. [DOI:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

### **Coverage Estimation**
- **CoverM**: https://github.com/wwood/CoverM

### **Metagenomic Virus Discovery**
- **geNomad**: https://github.com/apcamargo/genomad

## **Nextflow & Pipeline Framework**
- **Nextflow**: Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). "Nextflow enables reproducible computational workflows." *Nature Biotechnology*, 35(4), 316-319. [DOI:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## **Additional References**
For further reading, please refer to:
- **Nextflow Documentation**: [https://www.nextflow.io/docs/latest/](https://www.nextflow.io/docs/latest/)
- **Bioconda**: [https://bioconda.github.io/](https://bioconda.github.io/)

---

## **Acknowledgments**
This pipeline makes extensive use of publicly available tools, frameworks, and datasets. We thank the respective authors for making these tools and resources freely available to the scientific community.
