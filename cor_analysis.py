
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


def scan_tRNA(each_prefix_genome, each_genome_fastafile, output_path, each_organism = None):

    '''
    Scan tRNA in your organism of interest via 'tRNAscan-SE'.

    :param each_prefix_genome: prefix of genome fasta file passed from function getPrefix
    :param each_genome_fastafile: a fasta file of your organism of interest
    :param output_path: path/to/tRNAscan_count.txt
    :param organism: specify your organsim: bacteria, eukaryotic, archaeal and mitochondrial, otherwise the default is eukaryotic
    :return: a tRNAscan result txt file, with the tRNA copy number counted
    '''

    tRNAscan_out = join(output_path,each_prefix_genome+'_tRNAscan_count.txt')

    if each_organism == 'bacteria':
        mode = '-B'
    elif each_organism == 'eukaryotic':
        mode = '-E'
    elif each_organism == 'archaeal':
        mode = '-A'
    elif each_organism == 'mitochondrial':
        mode = '-M'
    elif each_organism is None:
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
    :return: df_tRNA_data, df_tRNAcount
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

    return tRNA_data_rc



def mergedf(output_path, each_prefix_cds):


    codon_count = pd.read_csv(join(output_path, each_prefix_cds + '_countBias_result.txt'), sep='\t')

    df1 = pd.DataFrame(codon_count, columns=['Codon', 'Bias'])
    df2 = pd.DataFrame(tRNA_count, columns=['Codon', 'tRNA_copy_N'])
    df1_2 = pd.merge(df1, df2, how='inner', on='Codon')
    df3 = pd.DataFrame(tRNA_data_rc, columns= ['Codon','fpkm'])
    df_1_2_3 = pd.merge(df1_2,df3, how='inner',on='Codon')
    # df_1_2_3 : codon,Bias,tRNA_copy_N,'fpkm'
    Df
    return df_1_2_3


def cor_analysis(df_1_2_3):



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

