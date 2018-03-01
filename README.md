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
  analysis on codon bias of each amino acid and tRNA copy numbers in your genomes.
- InputFile:CDS file in fasta format and genome file in fasta format. It support multiple input files.
- outputFile:your_prefix_count_result.txt; correlation_result.txt; correlation analysis plot

##
Workflow:
###
- If you start with raw reads, be sure to process it and assembly your reads into assemblies in fasta format using 
  some proper pipeline. Otherwise you can jump to next step directly.   
- Prepare your_organism_cds.fasta files using "transdecoder" or other similar CDS estimating packages.
- Run the codonCount.py script, see the usage manual below.
- Collect your codon count result downstream analysis.
- Searching for tRNA genes in genomic sequence.
- Correlation analysis on codon bias of each anmino acid and tRNA copy numbers.
- Save result and plot the correlation.
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
- software: tRNAscan-SE install via [here](http://lowelab.ucsc.edu/tRNAscan-SE/) and put into your path
##
Usage
```
Usage: codonCount.py 


arguments:
   -h, --help            show this help message and exit
  '-i', "--infile", help = "path/to/your/cds.fasta file, you must feed at least one cds.fasta file or multiple cds.fasta files, speperated by comma without space"
  '-o', "--outputDir", help="path/to/your/output_directory"
  '-p', "--prefix", help="a few characters of prefix to tag your output result files, e.g. sample1_; multiple prefixes should be passed correspondingly with your multiple input files,spereted by comma without space"

Usage: cor_tRNA_codon.py

arguments:
  '-i', "--genome_fastafile", help = "path/to/your/genome.fasta file, you must feed at least one genome.fasta file or multiple genome.fasta files, speperated by comma"
  '-c', "--codoncount_file", help = "path/to/your/codoncount_file, the result file of step 1. you must feed at least one genome.fasta file or multiple genome.fasta files, speperated by comma"
  '-o', "--outputDir", help = "path/to/your/output_directory"


```
