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

### * 1.Preprocess
  The presence of poor quality or technical sequences such as adapters in the sequencing data can easily result in suboptimal downstream analyses. There are many useful read preprocessing tools to perform the quality control (FastQC, Trimmomatic). Here we choose Trimmomatic to clean our sequencing datasets, for example:
  
    java -jar xx/software/trimmomatic-0.xx.jar PE -threads 30 -phred64 xx/R1.fq xx/R2.fq R1_paired_trimmed.fq R1_unpaired_trimmed.fq R2_paired_trimmed.fq R2_unpaired_trimmed.fq ILLUMINACLIP:xx/Trimmomatic-0.xx/adapters/TruSeqxxx.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20

<br>

### * 2.Assembly
  There are two types of sequence assembly, that is overlap/layout/consensus approach and de Bruijn graph approach. Because of the characteristics of next-generation sequencing platform and the incompletion of reference genome, de nove assembly is appropriate for metagenomics assembly. According to our experience, we use SPAdes Genomes Assembler to assembly our data. Executinng SPAdes as the following command:
  
    spades.py -o output -k 21,33,55 --pe1-1 xx/R1_paired_trimmed.fq --pe1-2 xx/R2_paired_trimmed.fq --pe1-s xx/R1R2_unpaired_trimmed.fq --mp1-1 xx/MP_R1_paired_trimmed.fq --mp1-2 xx/MP_R2_paired_trimmed.fq --mp1-s xx/MP_R1R2_unpaired_trimmed.fq --pacbio xx/filtered_subreads.fastq --careful --cov-cutoff 3 --disable-gzip-output -t 20 -m 500

<br>

### * 3.Binning
  In order to bin metagenomics sequence, we have to collect some fetures about the assembled fragments. The output of some assembler, such as SPAdes, Velvet and AByss, provide the coverage value of each cotigs in the header. If there isn’t coverage information in the assembly data, we can map reads to contigs using BWA or Bowtie to get coverage information instead. For each contigs, GC content, length and coverage are extracted from assembly file using my perl script. The script outputs the file contigs.info including 4 column, that is contig name, contig length, GC content and coverage. Additionly, we utilize the essential gene sets to draw back some assembly fragements from the unassigned sequences. The process of obtaining the information of essential gene contained in the assembly was doned by the same perl script.

    calc_contiginfo_essentialgene_v3.pl -r scaffolds.fasta -1 R1_paired_trimmed.fq -2 R2_paired_trimmed.fq -m 500 -I 0 -X 600 -p phred64 -t 30
Binning the metagenomic assembly using selected features

    Meta-binning.R
<br>

### * 4.Taxonomy
  We implemented taxonomic assignments of the genome bins using TAXAassign with some modiﬁed codes for efﬁciency and accuracy. We used deduced amino acid sequence information through DIAMOND BLASTP searches, instead of nucleotide sequences through BLASTN searches, to produce a protein sequence alignment against the NCBI non-redundant (nr) protein database. 
    
    TAXAassign_prot.sh -c 30 -r 100 -t 60 -m 60 -q 50 -a "60,65,70,80,95,95" -f All_bins.faa

<br>

### * 5.Bins Evaluation
  Since genomes reconstructed from metagenomic data usually vary substantially in their qualities, we proposed a set of quality criteria with quantitative thresholds to evaluate the quality of these genomes for subsequent analyses.

![image](https://github.com/qi-lee/Meta-Microbiome/blob/master/5_Bins_Evaluation/quality.JPG)
Other related software: CheckM.    
<br>

### * 6.Phylogenetic analysis
To infer phylogenetic relationships among bacteria, a whole-genome based and alignment-free Composition Vector Tree (CVTree) method was applied to the comparison and clustering of the 68 genomes we extracted from our assemblies.

    cvtree -i species.list -p data -o CVTree_k6.txt -k 6
    neighbor
Other related software: PHYLIP, MEGAN.
<br>

### * 7.Functional analysis

Functional characterization and annotation of protein encoding genes were performed by MOCAT2. To further compare the functional potential of each group, the predicted ORFs were analyzed using the GhostKOALA service on the KEGG website. When the results were returned through your email, open the links and then compare your interesting pathway.
    
    http://www.kegg.jp/ghostkoala/
Other related software: MEGAN.

<br>

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
