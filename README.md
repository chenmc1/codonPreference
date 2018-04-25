# CodonPreference User Manual 
####
Authors: Chenming Cui and Aureliano Bombarely
##
##
Citation: 
- Lowe, T.M. & Eddy, S.R. (1997) tRNAscan-SE: A program for improved detection of transfer RNA genes in genomic sequence. Nucl. Acids Res. 25: 955-964.
## 
Introduction:
###
- CodonPreference is designed to analysis codon bias/preference in your organisms of interest. Specially, it performs correlation
  analysis on codon bias of degenerated codons, tRNA gene copy numbers and its expression level in your genomes.
  Codon bias was defined as: (each condon count in the degenerate codons - average codon count of the degenerate codons)/average codon count of the degenerate codons.
- InputFile: CDS file in fasta format, genome file in fasta format, reads count matrix file with gene id and reads count columns header as "Sequence_Name" and "Count". It support multiple input files.
  Please specify your organsim from these categories: bacteria, eukaryotic, archaeal and mitochondrial, otherwise the default is eukaryotic.
  Reads count matrix file can be generated from RNASeq data using this script: https://github.com/chenmc1/illuminaRaw2CountMatrix-DGE
- OutputFile: your_prefix_count_result.txt; correlation_result.txt; correlation analysis plot.

##
Workflow:
###
- If you start with raw reads, be sure to process it and assembly your reads into assemblies in fasta format using 
  some proper pipeline. Otherwise you can jump to next step directly.   
- Prepare your_organism_cds.fasta files using "transdecoder" or other similar CDS estimating packages.
- Counted frequency of codons and its bias in the cds fasta file.
- Searched for tRNA genes in genomic sequence via tRNAscan-SE.
- Calculated the expression levels ("fpkm") of tRNA genes in your genome.  
- Correlation analysis on codon bias of each anmino acid, tRNA gene copy numbers and tRNA gene expression levels.
- Saved result and plot the correlation.
##
Dependencies and Requirements:
###
- python 2.x/3.x
- python modules (modules can be installed via [pip](https://pip.pypa.io/en/stable/installing/): pip install <module_name>)
  * BioPython
  * pandas
  * seaborn
  * matplotlib
  * numpy
  * csv
- software: tRNAscan-SE install via [here](http://lowelab.ucsc.edu/tRNAscan-SE/) and put it into your path
##
Usage


```
usage: wrapper.py [-h]
                  [-organism {bacteria,eukaryotic,archaeal,mitochondrial}]
                  genome_fastafile cds_fastafile reads_countfile outputDir

corelation analysis on condon bias, tRNA gene copy number and expression level
in the genome of your interest.

positional arguments:
  genome_fastafile      path/to/your/genome.fasta file, you must feed at least
                        one genome.fasta file or multiple genome.fasta files,
                        speperated by comma
  cds_fastafile         path/to/your/cds.fasta file', you must feed at least
                        one cds.fasta file or multiple genome.fasta files,
                        speperated by comma
  reads_countfile       path/to/your/reads_countfile, you must feed at least
                        one reads_countfile or multiple reads_count files,
                        speperated by comma
  outputDir             path/to/your/output_directory

optional arguments:
  -h, --help            show this help message and exit
  -organism {bacteria,eukaryotic,archaeal,mitochondrial}
                        Please specify the organism that you are searching
                        tRNA in these four organsim types
                        'bacteria','eukaryotic','archaeal','mitochondrial'.The
                        default organsim is eukaryotic'



```

## 
About
###
[Contact me!](mailto:chenmc1@vt.edu)



