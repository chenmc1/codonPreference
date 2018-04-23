#!/usr/bin/env python

__author__ = 'chenmingcui'


import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
import numpy as np
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
from __future__ import division
from config import check_dependencies,checkHeader_readscountfile,check_CDSfasta
import config
from codonCount import *
from cor_analysis import *


#### test dir /Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/test0318

def get_args():
    parser = argparse.ArgumentParser(
        description="corelation analysis on condon bias, tRNA gene copy number and expression level in the genome of your interest.")

    parser.add_argument("genome_fastafile",
                        help="path/to/your/genome.fasta file, you must feed at least one genome.fasta file or multiple genome.fasta files, speperated by comma")

    parser.add_argument("cds_fastafile",
                        help="path/to/your/cds.fasta file', you must feed at least one cds.fasta file or multiple genome.fasta files, speperated by comma")

    parser.add_argument("reads_countfile",
                        help="path/to/your/reads_countfile, you must feed at least one reads_countfile or multiple reads_count files, speperated by comma")

    parser.add_argument("outputDir", help="path/to/your/output_directory")

    parser.add_argument('-organism', default='eukaryotic',
                        choices=['bacteria', 'eukaryotic', 'archaeal', 'mitochondrial'],
                        help="Please specify the organism that you are searching tRNA in these four organsim types 'bacteria','eukaryotic','archaeal','mitochondrial'."
                             "The default organsim is eukaryotic'")  # optional argument

    args = parser.parse_args()

    return args



def Main():

    # item 1: pass arguments
    args = get_args()

    genome_fastafiles = args.genome_fastafile.split(",")
    cds_fastafiles = args.cds_fastafile.split(",")
    reads_countfiles = args.reads_countfile.split(",")
    organisms = args.organism.split(",")
    output_path = args.outputDir

    for each_genome_fastafile in genome_fastafiles:
        for each_cdsfile in cds_fastafiles:
            for each_readscountfile in reads_countfiles:
                for each_organism in organism:
                    Main(each_genome_fastafile, each_cdsfile, each_readscountfile, output_path, each_organism)


    # item 2: pass dictionaries
    dict1 = config.D1
    dict2 = config.D2
    dict3 = config.D3


    # item 3: check if tRNAscan-SE installed
    if check_dependencies('tRNAscan-SE') is False:
        print('tRNAscan-SE was not installed, please install and put it into your path; http://lowelab.ucsc.edu/tRNAscan-SE')
        sys.exit(2)
    else:
        print('tRNAscan-SE was installed successfully and ready for the following analysis.')


    # item 4: QC on your.cds.fasta file
    if check_CDSfasta(each_cdsfile) is False:
        print("your cds.fasta file didn't pass the quality check, make sure you generate cds.fasta file from transDecoder. https://github.com/TransDecoder/TransDecoder/wiki")
        sys.exit(3)
    else:
        print("Input-cds-file quality pre-check passed!")


    # item 5: check if reads count header named properly
    if checkHeader_readscountfile(each_readscountfile) is False:
        print(
            "File header Error: unable to find 'Sequence_Name' and 'Count' columns in your reads count file or you didn't name these two columns exactly as Sequence_Name and Count.")
        sys.exit(4)
    else:
        print('reads count file header passed check!')


    # item 6: get prefix from the files that user provides
    each_prefix_cds, each_prefix_genome = getPrefix(each_cdsfile,each_genome_fastafile)


    # item 7: read each_cdsfile
    seqs = readCDSfasta(each_cdsfile)


    # item 8: count codon on each of the transcripts(seqs)
    dict2 = countCodon(seqs)


    # item 9: Write codon count result into a table.txt file
    Codoncount2table(output_path, each_prefix_cds)


    # item 10: Scan tRNA in your organism of interest via 'tRNAscan-SE'
    scan_tRNA(each_prefix_genome, each_genome_fastafile, output_path, each_organism=None)


    # item 11: Process tRNAscan-SE results
    df_tRNA_data, df_tRNAcount = process_tRNAscan(output_path,each_prefix_genome)


    # item 12: Read and process reads count file
    df_rc, totalRC = process_readsCount(each_readscountfile)


    # item 13: Merge the two data frame to get the tRNA gene reads count from each reads count file
    tRNA_data_rc = mergeGet64count(df_tRNA_data, df_rc)


    # item 14: Calculate fpkm and append fpkm to tRNA_data_rc
    tRNA_data_rc = tRNAgene_fpkm(tRNA_data_rc, totalRC, output_path, each_prefix_genome)


    # item 15: Merge to get a data frame, df_1_2_3 : codon,Bias,tRNA_copy_N,'fpkm'
    df_1_2_3 = mergedf(output_path, each_prefix_cds)


    # item 16: Correlation analysis between on df_1_2_3 and write results into a file
    cor_analysis(df_1_2_3)


    # item 17: plot the correlationship curve and save it to your output dir
    plot(df_1_2_3)


if __name__ =='__main__':
    if len(sys.argv) ==1 :
         print ("please provide 4 necessary files to run the program!")
         print ("please refer to the user manual or -h for help")

    Main()























