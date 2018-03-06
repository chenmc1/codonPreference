#!/usr/bin/python
__author__ = 'chenmingcui'

import os,sys
import argparse
from os.path import join
import pandas as pd
from Bio.Seq import Seq
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#    cor      https://stackoverflow.com/questions/3949226/calculating-pearson-correlation-and-significance-in-python

def cor_tRNA_codon(each_genome_fastafile,output_path,each_codoncount_file):
    """"
    This function will perform the correlation analysis between tRNA copy numbers and count count number
    step 1 : count tRNA from genome
    step 2 : correlation analysis
    """
    dict_3 = {'TTT': 0, 'TTC': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0,
              'TAT': 0, 'TAC': 0, 'TGT': 0, 'TGC': 0, 'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0,
              'CTG': 0, 'TTA': 0, 'TTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'ATT': 0,
              'ATC': 0, 'ATA': 0, 'ATG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0,
              'AGG': 0, 'TTT': 0, 'TTC': 0, 'TTT': 0, 'TTC': 0, 'ACT': 0, 'ACC': 0, 'ACA': 0,
              'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 'TTT': 0, 'TTC': 0, 'CAA': 0, 'CAG': 0,
              'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0,
              'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,
              'CAT': 0, 'CAC': 0, 'TAA': 0, 'TAG': 0, 'TGA': 0, 'TGG': 0, 'ACG': 0, "NNN": 0}

    ## count tRNA from genome using software tRNAscan
    each_prefix,file_ext =os.path.splitext(each_genome_fastafile)
    cmd = "tRNAscan-SE "+each_genome_fastafile+" -o "+join(output_path,each_prefix+'_tRNA_count.txt')
    print (cmd)
    os.system(cmd)

    ## process the tRNA scan result, read it as dataframe in pandas
    #tRNA_data = pd.read_csv('Cannuum_CA00562_Trinity_tRNA_count.txt', sep='\t', skiprows=(0, 1), header=0)
    tRNA_data = pd.read_csv(join(output_path, each_prefix + '_tRNAscan_result.txt'), sep='\t', skiprows=(0,1), header=0)
    anticodon_df = tRNA_data.iloc[:, 5]  # slice the anti_codon column #5
    anticodon_list = list(anticodon_df)  # anticodons are in the list now
    tRNA_data['length'] = tRNA_data.iloc[:, 2].sub(tRNA_data.iloc[:, 3], axis=0).abs()# getting gene length by substracting
                                                                                     # col2 and col3 based on column index0


    ## complement and reverse the anticodon for matching the codon dictionary
    my_codon_list = []
    for anticodon in anticodon_list:
        my_dna = Seq(anticodon)
        my_codon = str(my_dna.reverse_complement())  # append only work for str
        my_codon_list.append(my_codon)

    ## count the "codon" and update the dict_3 with the counting
    for my_codon in my_codon_list:
        dict_3[my_codon] +=1

    ## write the tRNA copy number into a csv file
    fieldnames = ['Codon', 'tRNA_copy_N'] # the result of tRNA count numbers in your genome
    with open(join(output_path, each_prefix + '_anticodonCount.csv'), 'w') as f:  # anticodonCount is the tRNA numbers
        writer = csv.writer(f)
        writer.writerow(fieldnames)
        for k, v in dict_3.items():
           writer.writerow([k,v])

    ## correlation analysis between codon bias and tRNA copy numbers
    codon_count = pd.read_csv(each_codoncount_file, sep='\t')
    tRNA_count = pd.read_csv(join(output_path, each_prefix + '_anticodonCount.csv'))
    #tRNA_count = pd.read_csv('anticodonCount.csv')
    df1 = pd.DataFrame(codon_count, columns=['Codon', 'Bias'])
    df2 = pd.DataFrame(tRNA_count, columns=['Codon', 'tRNA_copy_N'])
    merged_inner = pd.merge(df1, df2, how='inner', on='Codon')
    merged_inner.to_csv(join(output_path,each_prefix + 'merged_codon_tRNA.csv'))
    cor = merged_inner['Bias'].corr(merged_inner['tRNA_copy_N'])

    ## save correlation to file and plot the correlation
    final = open(join(output_path, each_prefix + '_correlation.txt'), 'w')
    final.write(cor)
    final.close ()
    sns.regplot(x="Bias", y="tRNA_copy_N", data=merged_inner)
    #plt.show()  # plot in the pycharm
    plt.savefig(join(output_path, each_prefix + '_correlation.plot')


def Main():
    parser = argparse.ArgumentParser(description= "corelation analysis between condon bias of each anmino acids with tRNA copy numbers in the genome")

    parser.add_argument('-i', "--genome_fastafile", help = "path/to/your/genome.fasta file, you must feed at least one genome.fasta file or multiple genome.fasta files, speperated by comma")
    parser.add_argument('-c', "--codoncount_file", help = "path/to/your/codoncount_file, the result file of step 1. you must feed at least one genome.fasta file or multiple genome.fasta files, speperated by comma")
    parser.add_argument('-o', "--outputDir", help = "path/to/your/output_directory")

    args = parser.parse_args()

    genome_fastafiles = args.genome_fastafile.split(",")
    codoncount_files = args.codoncount_file.split(",")
    output_path = args.outputDir

    for each_genome_fastafile in genome_fastafiles:
           for each_codoncount_file in codoncount_files:
                 cor_tRNA_codon(each_genome_fastafile, output_path,each_codoncount_file)      ## function_name(argv.pass_variable)



if __name__ =='__main__':
    Main()

## test code
# codon_count = pd.read_csv('SAM1__count_result.txt', sep = '\t')
# tRNA_count = pd.read_csv('anticodonCount.csv')
# df1 = pd.DataFrame(codon_count, columns=['Codon', 'Bias'])
# df2 = pd.DataFrame(tRNA_count, columns=['Codon', 'tRNA_copy_N'])
# merged_inner = pd.merge(df1, df2, how='inner', on='Codon')
# cor = merged_inner['Bias'].corr(merged_inner['tRNA_copy_N'])
#

