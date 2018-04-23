#!/usr/bin/python
from __future__ import division
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

#### test dir /Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference/test0318

dict_1 = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
          "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
          "TAT": "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP",
          "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
          "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
          "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
          "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
          "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
          "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
          "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

dict_2 = {
    'F': {'TTT': 0, 'TTC': 0}, 'S': {'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0},
    'Y': {'TAT': 0, 'TAC': 0}, 'C': {'TGT': 0, 'TGC': 0}, 'W': {'TGG': 0},
    'L': {'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'TTA': 0, 'TTG': 0},
    'P': {'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0},
    'I': {'ATT': 0, 'ATC': 0, 'ATA': 0}, 'M': {'ATG': 0},
    'R': {'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0},
    'F': {'TTT': 0, 'TTC': 0}, 'F': {'TTT': 0, 'TTC': 0}, 'T': {'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0},
    'N': {'AAT': 0, 'AAC': 0}, 'K': {'AAA': 0, 'AAG': 0}, 'F': {'TTT': 0, 'TTC': 0}, 'Q': {'CAA': 0, 'CAG': 0},
    'V': {'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0}, 'A': {'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0},
    'G': {'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}, 'D': {'GAT': 0, 'GAC': 0}, 'E': {'GAA': 0, 'GAG': 0},
    'H': {'CAT': 0, 'CAC': 0},
    'STOP': {'TAA': 0, 'TAG': 0, 'TGA': 0, 'TGG': 0}
}

dict_3 = {'TTT': 0, 'TTC': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0,
          'TAT': 0, 'TAC': 0, 'TGT': 0, 'TGC': 0, 'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0,
          'CTG': 0, 'TTA': 0, 'TTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ATT': 0,
          'ATC': 0, 'ATA': 0, 'ATG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0,
          'AGG': 0, 'TTT': 0, 'TTC': 0, 'TTT': 0, 'TTC': 0, 'ACT': 0, 'ACC': 0, 'ACA': 0,
          'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'TTT': 0, 'TTC': 0, 'CAA': 0, 'CAG': 0,
          'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
          'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
          'CAT': 0, 'CAC': 0, 'TAA': 0, 'TAG': 0, 'TGA': 0, 'TGG': 0, 'ACG': 0, "NNN": 0}


each_cdsfile ='Cannuum_CA00562_Trinity.cds.fa'
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
    check module funtion by help(readCDSfasta)
    1. check file existence
    2. read each_cdsfile, it should be the cds fasta file that generated from 'transdecoder'
    3. check quality of the cds fasta file

    :param each_cdsfile: input cds fasta file
    :return: pass or warning
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

    ## input file quality check
    for k, v in seqs.items():
        # print(k,v)
        if len(v) % 3 != 0:
            print(
                "Warning: You have at least one transcript(%s) that is not in reading frame, be sure to estimate cds using TransDecoder" %k)
            break
        elif v[-3:].upper() not in ["TAA", "TAG", "TGA"]:
            print(
                "Warning: You have at least one transcript: %s that is not terminated by stop codon, be sure to estimate cds using TransDecoder" %k)
            break
        elif Seq(v).translate().find("*") < len(v) / 3 - 1:
            print(
                "Warning: You have at least one transcript: %s that has a premature stop codon" %k)
            break

        else:
            print("Input-cds-file quality pre-check passed!")
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

    :param output_path: define your output path
    :param each_prefix_cds: pass prefix that parsed from getPrefix
    :return: a table with four column that collects the codon usage frequency (count) and percentage (Bias) for each amino acids.
    '''
    # this_dict1 = dict_1
    # this_dict2 = dict_2
    f = open(join(output_path, each_prefix_cds + '_countBias_result.txt'), 'w')
    f.write(("AminoAcid\tCodon\tCount\tBias\n"))
    codon_count_pair = list(dict2.values())  # values are in a list [{'GAT': 78, 'GAC': 18}, {'GTG': 38, 'GTA': 20, 'GTT': 50, 'GTC': 24}]
    for n in range(0, len(codon_count_pair)):
        codons = list(codon_count_pair[n].keys())
        for m in range(0, len(codons)):
            count = codon_count_pair[n][codons[m]]
            bias = float(count) / sum(codon_count_pair[n].values())
            f.write(dict_1[codons[m]] + '\t' + codons[m] + '\t' + str(count) + '\t' + '{0:.2f}'.format(bias) + '\n')
    f.close()

    print('Counting of codon finished, check your output file.')



def scan_tRNA(each_prefix_genome, each_genome_fastafile, output_path, organism = None):

    '''
    Scan tRNA in your organism of interest via 'tRNAscan-SE'.

    :param each_prefix_genome: prefix of genome fasta file passed from function getPrefix
    :param each_genome_fastafile: a fasta file of your organism of interest
    :param output_path: path/to/tRNAscan_count.txt
    :param organism: specify your organsim: bacteria, eukaryotic, archaeal and mitochondrial, otherwise the default is eukaryotic
    :return: a tRNAscan result txt file, with the tRNA copy number counted
    '''

    tRNAscan_out = join(output_path,each_prefix_genome+'_tRNAscan_count.txt')

    if organism == 'bacteria':
        mode = '-B'
    elif organism == 'eukaryotic':
        mode = '-E'
    elif organism == 'archaeal':
        mode = '-A'
    elif organism == 'mitochondrial':
        mode = '-M'
    elif organism is None:
        mode = '-E'
    else:
        print('scan_tRNAError: Please specify the organism that you are searching tRNA in these four organsim types"bacteria","eukaryotic","archaeal","mitochondrial".The default organsim is eukaryotic')

        sys.exit(0)

    subprocess.call(['tRNAscan-SE', each_genome_fastafile, mode, '-o', tRNAscan_out])


def process_tRNAscan(output_path,each_prefix_genome):
    # tRNA_data = pd.read_csv('Cannuum_CA00562_Trinity_tRNA_count.txt', sep='\t', skiprows=(0, 1), header=0)

    '''
    1. Read tRNAscan result file in a dataframe: tRNA_data
    2. complementary and reverse the anticodon in column 5, count the number
    :param output_path: output_path
    :param each_prefix_genome: the prefix
    :return:
    '''
    ## process the tRNA scan result, read it as dataframe in pandas
    tRNAscan_out = join(output_path, each_prefix_genome + '_tRNAscan_count.txt')
    tRNA_data = pd.read_csv(tRNAscan_out, sep='\t', skiprows=(0, 1),header=0)

    ## complement and reverse the anticodon for matching the codon dictionary
    anticodon_df = tRNA_data.iloc[:, 5]  # slice the anti_codon column #5
    anticodon_list = list(anticodon_df)  # anticodons are in the list now
    my_codon_list = []
    for anticodon in anticodon_list:
        my_dna = Seq(anticodon)
        my_codon = str(my_dna.reverse_complement())  # append only work for str
        my_codon_list.append(my_codon)

    ## append codon list to tRNA_data
    se = pd.Series(my_codon_list)
    tRNA_data['codon'] = se.values
    df_tRNA_data = pd.DataFrame(tRNA_data)

    ## count the "codon" and update the dict_3 with the counting
    for my_codon1 in my_codon_list:
        dict_3[my_codon1] +=1

    ## write Codon and tRNA_copy_number into a dataframe df_tRNAcount
    df_tRNAcount = pd.DataFrame()
    df_tRNAcount['Codon'] = dict_3.keys()
    df_tRNAcount['tRNA_copy_N'] = dict_3.values()


    return df_tRNA_data, df_tRNAcount
    #df_tRNA_data,df_tRNAcount = process_tRNAscan()

    # !! tRNA_data is a dataframe has sequence name, length and codon columns etc. now


def process_readsCount(each_readscountfile):

    '''
    reading a reads count file in csv format, the column with gene id
    :param each_readscountfile: reads count file that generated from RNASeq data
    :return: dataframe(df_rc) that has two columns "Sequence_Name" and "Count"
    '''
    ## read the read_count file as dataframe in pandas
    readsCount_df = pd.read_csv(each_readscountfile)
    df_rc = pd.DataFrame(readsCount_df, columns=['Sequence_Name', 'Count'])
    totalRC = df_rc.sum(axis=1)  ## total reads count in the count file, (maybe sum of 30000 genes' reads)

    return df_rc,totalRC

  #  df_rc, totalRC = process_readsCount()



def mergeGet64count(df_tRNA_data, df_rc):
    '''
    merge the two data frame to get the tRNA gene reads count from each reads count file.
    :return: tRNA_data_rc has columns: Sequence name, length, codon, count.   64 genes rather than all genes in count files
    '''
    ### merge the reads count dataframe and the tRNA_data dataframe based on the sequence name
    tRNA_data_rc = pd.merge(df_tRNA_data, df_rc, how='inner', on='Sequence Name')
    if tRNA_data_rc.empty:
        print('Dateframe merge Error! Dataframe: tRNA_data_rc is empty.')
        print('Making sure that 'Sequence Name' in reads count file and tRNA-scan result file is same sort of gene IDs')
        sys.exit(1)

    return tRNA_data_rc




def tRNAgene_fpkm(tRNA_data_rc, totalRC, output_path,each_prefix_genome):

    tRNAscan_out = join(output_path, each_prefix_genome + '_tRNAscan_count.txt')
    tRNA_data = pd.read_csv(tRNAscan_out, sep='\t', skiprows=(0, 1), header=0)
    tRNA_data['length'] = tRNA_data.iloc[:, 2].sub(tRNA_data.iloc[:, 3], axis=0).abs()  # add a column: gene length; getting gene length by substracting

    tRNA_data_rc['fpkm'] = tRNA_data_rc.iloc[:, 10] / tRNA_data_rc.iloc[:, 12] / totalRC  # col 12 is the gene length





def mergedf(output_path, each_prefix_cds):


    codon_count = pd.read_csv(join(output_path, each_prefix_cds + '_countBias_result.txt'), sep='\t')

    df1 = pd.DataFrame(codon_count, columns=['Codon', 'Bias'])
    df2 = pd.DataFrame(tRNA_count, columns=['Codon', 'tRNA_copy_N'])
    df1_2 = pd.merge(df1, df2, how='inner', on='Codon')
    df3 = pd.DataFrame(tRNA_data_rc, columns= ['Codon','fpkm'])
    df_1_2_3 = pd.merge(df1_2,df3, how='inner',on='Codon')
    # df_1_2_3 : codon,Bias,tRNA_copy_N,'fpkm'

    return df_1_2_3


def cor_analysis(Df):

    Df = df_1_2_3

    df_1_2_3.to_csv(join(output_path, each_prefix_cds + 'stats4cor.csv'))
    cor_Bias_N = df_1_2_3['Bias'].corr(df_1_2_3['tRNA_copy_N'])
    cor_Bias_exp = df_1_2_3['Bias'].corr(df_1_2_3['fpkm'])
    with open(join(output_path, each_prefix_cds + 'stat4cor.csv'), 'a') as f:
        f.write('correlation coefficient of codon bias and tRNA copy number is:%d\n' % (cor_Bias_N),
                'correlation coefficient of codon bias and tRNA expression level is:%d' % (cor_Bias_exp))


def plot(df_1_2_3):

    sns.regplot(x="Bias", y="tRNA_copy_N", data=df_1_2_3,
                scatter_kws={'color': 'black'}, line_kws={'color': 'red'})
    # plt.show()  # plot in the pycharm
    plt.savefig(join(output_path, each_prefix_cds + '_corB_N.plot'))
    plt.clf()
    sns.regplot(x="Bias", y="fpkm", data=df_1_2_3,
                scatter_kws={'color': 'black'}, line_kws={'color': 'blue'})
    # plt.show()  # plot in the pycharm
    plt.savefig(join(output_path, each_prefix_cds + '_corB_exp.plot'))





######
##### test









