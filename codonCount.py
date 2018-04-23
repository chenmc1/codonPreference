

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


#each_cdsfile ='Cannuum_CA00562_Trinity.cds.fa'
def getPrefix(each_cdsfile,each_genome_fastafile):

    '''
    Program will automatically get prefix from user's input file name to tag all the downstream resulting files.
    Please name your cds and genome files properly with the rules of: sample_prefix_cds.fasta or sample_prefix_genome.fasta

    :param each_cdsfile: input cds.fasta file
    :param each_genome_fastafile: input genome fasta file
    :return: the prefix - periods that besides last period of the file name will be the prefix
             for example: prefix of Cannuum_CA00562_Trinity.cds.fa is Cannuum_CA00562_Trinity.cds
    '''

    each_prefix_cds, file_ext1 = os.path.splitext(each_cdsfile)
    each_prefix_genome, file_ext2 = os.path.splitext(each_genome_fastafile)
    return each_prefix_cds,each_prefix_genome

def readCDSfasta(each_cdsfile):

    '''

    read each_cdsfile, it should be the cds fasta file that generated from 'transdecoder'


    :param each_cdsfile: input cds fasta file
    :return: seqs
    '''

    try:
        f = open(each_cdsfile)
    except IOError:
        print('No such file: %s' % each_cdsfile)

    ## read a fasta file
    seqs = {}
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            words = line.split()
            name = words[0][1:]
            seqs[name] = ''
        else:
            seqs[name] = seqs[name] + line
    f.close()


    return seqs



def countCodon(seqs):

    '''
    count codon on each of the transcripts(seqs)

    :param seqs: seqs are read and collected from cds fasta file
    :return: the counted dict_2
    '''

    seq = list(seqs.values())
    for i in range(0, len(seq)):
        each_seq = seq[i]  # seq is a list of transcripts
        for s in range(0, len(each_seq), 3):  # fetch each of the 3 nt in your transcript
            codon = each_seq[s:s + 3].upper()
            aa = dict1[codon]
            codon_count_dict = dict2[aa]
            codon_count_dict[codon] += 1
    return dict2



def Codoncount2table(output_path, each_prefix_cds):

    '''
    Write codon count result into a table.txt file.

    codonBias was defined as: (each condon count in the degenerate codons - average codon count of the degenerate codons)/average codon count of the degenerate codons

    :param output_path: define your output path
    :param each_prefix_cds: pass prefix that parsed from getPrefix
    :return: a table with four column that collects the codon usage frequency (count) and percentage (Bias) for each amino acids.
    '''

    f = open(join(output_path, each_prefix_cds + '_countBias_result.txt'), 'w')
    f.write(("AminoAcid\tCodon\tCount\tBias\n"))
    codon_count_pair = list(dict2.values())  # values are in a list [{'GAT': 78, 'GAC': 18}, {'GTG': 38, 'GTA': 20, 'GTT': 50, 'GTC': 24}]
    for n in range(0, len(codon_count_pair)):
        codons = list(codon_count_pair[n].keys())
        # len(codons) : number of degenerate codons
        for m in range(0, len(codons)):
            count = codon_count_pair[n][codons[m]] #count on each codon
            average_count = sum(codon_count_pair[n].values())/len(codons) # average count for each set of degenerate codons
            bias0 = float(count)-(average_count) /average_count
            bias = abs(bias0)
            f.write(dict_1[codons[m]] + '\t' + codons[m] + '\t' + str(count) + '\t' + '{0:.2f}'.format(bias) + '\n')
    f.close()

    print('Counting of codon finished, check your output file.')


    ##