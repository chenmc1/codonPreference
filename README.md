# CodonPreference User Manual 
####
Authors: Chenming Cui and Aureliano Bombarely
## 
Introduction:
###
codonPreference is designed to analysis codon bias/preference in your organisms of interest.
I will keep it update for the downstream related analysis.

InputFile: CDS file in fasta format;
outputFile: your_prefix_count_result.txt
##
Workflow:
###
- If you start with raw reads, be sure to process it and assembly your reads into assemblies in fasta format using 
  some proper pipeline. Otherwise you can jump to next step directly.   
- Prepare your_organism_cds.fasta files using "transdecoder" or other similar CDS estimating packages.
- Run the codonCount.py script, see the usage manual below.
- Collect your codon count result.
##
Dependencies and requirements
###
- python 2.x/3.x
- BioPython v1.64 above
##
Usage
```
Usage: codonCount.py [-h] [--outputDir OUTPUTDIR] [--prefix PREFIX] infile

positional arguments:
  infile                path/to/your/cds.fasta file, you must feed at least
                        one cds.fasta file or multiple cds.fasta files,
                        speperated by comma

optional arguments:
  -h, --help            show this help message and exit
  --outputDir OUTPUTDIR
                        path/to/your/output_directory
  --prefix PREFIX       a few characters of prefix to tag your output result
                        files, e.g. gsf533_; multiple prefixes should be
                        passed correspondingly with your multiple input
                        files,spereted by comma

```
