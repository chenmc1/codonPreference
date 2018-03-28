
#!/usr/bin/python
from __future__ import division
__author__ = 'chenmingcui'


import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
import numpy as np


#filename = sys.argv[1:]
#outpath = sys.argv[2]
#def codonCount("/Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/sample1.cds.fasta"):
def codonCount(each_infile,output_path,each_prefix):
    """
    This function is to count codon useage in your transcriptome.
    The input file is the cds fasta file that generated from transdecoder.
    The ouput file is a summary of count for each codon in a tab delimited text file.

    :param filename:
    :return:
    """
    ## dictionary 1 : condons are translated to amino acids
    dict_1 = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
       "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    ## dictionary 2 : initialize codon counts, 0 by default
    dict_2 = {
          'F':{'TTT':0,'TTC':0},'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
          'Y':{'TAT':0,'TAC':0},'C':{'TGT':0,'TGC':0},'W':{'TGG':0},
          'L':{'CTT':0,'CTC':0,'CTA':0,'CTG':0,'TTA':0,'TTG':0},'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
          'I':{'ATT':0,'ATC':0,'ATA':0},'M':{'ATG':0},'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
          'F':{'TTT':0,'TTC':0},'F':{'TTT':0,'TTC':0},'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
          'N':{'AAT':0,'AAC':0},'K':{'AAA':0,'AAG':0},'F':{'TTT':0,'TTC':0},'Q':{'CAA':0,'CAG':0},
          'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
          'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0},'D':{'GAT':0,'GAC':0},'E':{'GAA':0,'GAG':0},'H':{'CAT':0,'CAC':0},
          'STOP':{'TAA':0,'TAG':0,'TGA':0,'TGG':0}
    }

    ## check if file exists

    try:
       f = open(each_infile)
    except IOError:
       print('No such file: %s' %each_infile)

    ## read a fasta file
    seqs={}
    for line in f:
        line=line.rstrip()
        if line.startswith('>'):
            words=line.split()
            name=words[0][1:]
            seqs[name]=''
        else:
            seqs[name]=seqs[name]+line
    f.close()

    ## input file quality check
    for k,v in seqs.items():
        #print(k,v)
        if len(v) % 3 != 0:
            print("Warning: You have at least one transcript(%s) that is not in reading frame, be sure to estimate cds using TransDecoder" %k)
            break
        elif v[-3:].upper() not in ["TAA","TAG","TGA"]:
            print("Warning: You have at least one transcript: %s that is not terminated by stop codon, be sure to estimate cds using TransDecoder" %k)
            break
        elif Seq(v).translate().find("*") < len(v)/3-1:
            print("Warning: You have at least one transcript: %s that has a premature stop codon" %k)
            break

        else:
            print("Input file quality pre-check passed!")

    ## now, seqs is a dictionary having all the seq_id and seq_sequence
    ## access to each transcript and count the codons
    seq =list(seqs.values())   # convert dictionary values into list otherwise cannot be indexed
    for i in range (0,len(seq)):
        each_seq = seq[i]      # seq is a list of transcripts
        for s in range (0,len(each_seq),3):  # fetch each of the 3 nt in your transcript
            codon = each_seq[s:s+3].upper()
            aa = dict_1[codon]
            codon_count_dict = dict_2[aa]
            codon_count_dict[codon] += 1

    #for k, v in dict_2.items():   # dict2 now has the count updated, print to check
    #    print (k,v)


    ## print the count into a table

    f = open(join(output_path, each_prefix+'_count_result.txt'),'w')
    f.write(("AminoAcid\tCodon\tCount\tBias\n"))
    codon_count_pair = list(dict_2.values()) # values are in a list   [{'GAT': 78, 'GAC': 18}, {'GTG': 38, 'GTA': 20, 'GTT': 50, 'GTC': 24}]
    for n in range (0,len(codon_count_pair)):
           codons =list(codon_count_pair[n].keys())
           for m in range (0,len(codons)):
               count = codon_count_pair[n][codons[m]]
               bias = float(count)/sum(codon_count_pair[n].values())
               f.write(dict_1[codons[m]] + '\t' + codons[m] + '\t' + str(count) + '\t' + '{0:.2f}'.format(bias) + '\n')
    f.close()

    print('Counting finished, check your output file.')


def Main():
    parser = argparse.ArgumentParser(description= "calculate codon bias of a genome")

    parser.add_argument('-i', "--infile", help = "path/to/your/cds.fasta file, you must feed at least one cds.fasta file or multiple cds.fasta files, speperated by comma without space")
    parser.add_argument('-o', "--outputDir", help="path/to/your/output_directory")
    parser.add_argument('-p', "--prefix", help="a few characters of prefix to tag your output result files, e.g. sample1_; multiple prefixes should be passed correspondingly with your multiple input files,spereted by comma without space")

    args = parser.parse_args()

    infiles = args.infile.split(",")
    output_path = args.outputDir
    prefixes = args.prefix.split(",")
    for each_infile in infiles:
        for each_prefix in prefixes:
            codonCount(each_infile,output_path,each_prefix)      ## function_name(argv.pass_variable)



if __name__ =='__main__':
    Main()


#commend to run it
#python codonCount.py -i /Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/sample1.cds.fasta,/Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/sample2.cds.fasta,/Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/sample3.cds.fasta -o /Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/test_dir/ -p SAM1_,SAM2_,SAM3_
#
#