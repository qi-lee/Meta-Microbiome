# Meta-Microbiome analysis pipeline

## This package is built around a collection of publicly available tools and personal scripts tied together for analyzing metagenomic microbiome datasets.

### The pipeline was developed by QiLi (liqi at ihb.ac.cn, Lab of Algal Genomics).<br> External software used by this pipeline are copyright respective authors.
<br>

### The pipeline can be broadly separated into seven main sections：
* 1.Preprocess
* 2.Assembly
* 3.Binning
* 4.Taxonomy
* 5.Bins_Evaluation
* 6.Phylogenetic_analysis
* 7.Functional_analysis
<br>

### 1.Preprocess
  The presence of poor quality or technical sequences such as adapters in the sequencing data can easily result in suboptimal downstream analyses. There are many useful read preprocessing tools to perform the quality control (FastQC, Trimmomatic). Here we choose Trimmomatic to clean our sequencing datasets, for example:
  
    java -jar xx/software/trimmomatic-0.xx.jar PE -threads 30 -phred64 xx/R1.fq xx/R2.fq R1_paired_trimmed.fq R1_unpaired_trimmed.fq R2_paired_trimmed.fq R2_unpaired_trimmed.fq ILLUMINACLIP:xx/Trimmomatic-0.xx/adapters/TruSeqxxx.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20
     


### Reference:
    Andrews S. FastQC: a quality control tool for high throughput sequence data[J]. 2010.
    Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing[J]. Journal of computational biology, 2012, 19(5): 455-477.
    Bolger A M, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data[J]. Bioinformatics, 2014, 30(15): 2114-2120.
    Buchfink B, Xie C, Huson D H. Fast and sensitive protein alignment using DIAMOND[J]. Nature methods, 2015, 12(1): 59-60.
    Huson D H, Beier S, Flade I, et al. MEGAN community edition-interactive exploration and analysis of large-scale microbiome sequencing data[J]. PLoS computational biology, 2016, 12(6): e1004957.
    Ijaz U, Quince C. TAXAassign v0. 4[J]. 2013.
    Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes[J]. Nucleic acids research, 2000, 28(1): 27-30.
    Kanehisa M, Sato Y, Morishima K. BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences[J]. Journal of molecular biology, 2016, 428(4): 726-731.
    Kultima J R, Coelho L P, Forslund K, et al. MOCAT2: a metagenomic assembly, annotation and profiling framework[J]. Bioinformatics, 2016, 32(16): 2520-2523.
    Parks D H, Imelfort M, Skennerton C T, et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes[J]. Genome research, 2015, 25(7): 1043-1055.
    Peng Y, Leung H C M, Yiu S M, et al. Meta-IDBA: a de Novo assembler for metagenomic data[J]. Bioinformatics, 2011, 27(13): i94-i101.
    Plotree D, Plotgram D. PHYLIP-phylogeny inference package (version 3.2)[J]. cladistics, 1989, 5(163): 6.
    Qi J, Luo H, Hao B. CVTree: a phylogenetic tree reconstruction tool based on whole genomes[J]. Nucleic acids research, 2004, 32(suppl_2): W45-W47.
    Quinlan A R, Hall I M. BEDTools: a flexible suite of utilities for comparing genomic features[J]. Bioinformatics, 2010, 26(6): 841-842.
    Xie M, Ren M, Yang C, et al. Metagenomic analysis reveals symbiotic relationship among bacteria in microcystis-dominated community[J]. Frontiers in microbiology, 2016, 7.
<br>

#### Questions or Comments, please contact: liqi at ihb.ac.cn.
