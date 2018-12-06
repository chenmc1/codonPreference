#!/usr/bin/env python3

import os
import sys
import pandas as pd
from os.path import join
import glob
import scipy.stats as stats

def file_numbers(out_path):

    n = len(glob.glob1(out_path, '*_combined_codonbias.txt'))
    return n


#out_path = '/Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova'
def process_combined_sample(out_path):
    """

    :param out_path:
    :return: generate dataframes ('_allsamples_bias.txt') for each population, each dataframe has all sample's codon bias
    (each column represent one sample)
    """

    # populationA_combined_codonbias.txt

    for file in os.listdir(out_path):

        if file.endswith('_combined_codonbias.txt'):
            file_path, file_name = os.path.split(file)
            prefix, middle, file_ext = file_name.split('_')
            df_codonbias = pd.DataFrame()
            df_combine = pd.read_csv(join(out_path, file), sep='\t')

            codon = pd.DataFrame(df_combine.iloc[:,[2]])
            df_codonbias = pd.concat([df_codonbias,codon],ignore_index=True,axis=1) # append column to an existing dataframe

            n = file_numbers(out_path)

            for i in range(0, n):
                each_col = pd.DataFrame(df_combine.iloc[:,[4 + i * 4]])
                df_codonbias = pd.concat([df_codonbias,each_col],ignore_index= True,axis=1)

            df_codonbias.to_csv(join(out_path, prefix + "_allsamples_bias.txt"),sep="\t" )


def CondonBias_df_anova(out_path):

    """

    "Differential codon bias analysis"

    This function was designed to perform variance analysis of codon bias across your treatments/stage of interest on your genome.
    Specially.  One-way and up to two-way analysis of variance (ANOVA) can be performed by calling this function.

    1. read sample._countBias_result.txt
    AminoAcid	Codon	Count	Bias
        A	     GCA	12	    41.38%
        A	     GCC	4	    13.79%
        A	     GCT	12	    41.38%
        A	     GCG	1	     3.45%
        C	     TGC	4	    33.33%
        C 	     TGT	8	    66.67%
        E	     GAG	16	    35.56%
        E	     GAA	29	    64.44%
        D	     GAT	39	    81.25%
        D	     GAC	9	    18.75%
        G	     GGT	22	    37.93%

    2. filter the df, generate a new df with the most biased codon for each AminioAcid:
        A ~ 41.38
        C ~ 66.67
        E ~ 64.44
        ..

    3. anova analysis： e.g. time course study

      or comparision in multiple populations

                               | population1  | population2 | population3  | population4   ...
                   ---------------------------------------------------------------------------
                   lysine      | sample1.2.3..|sample1.2.3..|sample1.2.3.. |sample1.2.3..
                   proline     | sample1.2.3..|sample1.2.3..|sample1.2.3.. |sample1.2.3..
                   ...
                   ...
    4. Or two-way anova analysis:    w/  reps:
                                     w/o reps:
                                ..

    5. Or 2k factorial design:       w/  reps:
                                     w/o reps:

    input: multiple files with path sample.bias_count_result.txt

    :return:
    """
    #1. read codon_bias file
    #codon_count = pd.read_csv(join(output_path, each_prefix_cds + '_countBias_result.txt'), sep='\t')

    ###### $$$$$$$$ optional $$$$$$$$$$

    #2. filter the df
 #   df_aaBias = pd.DataFrame(codon_count, columns = ['AminoAcid','Bias'])
    #2.1 purge the larget biasd codon
 #   df_aaBias_max = df_aaBias.groupby(['AminoAcid'], sort=False)['Bias'].max()

    #       AminoAcid	Codon	Count	Bias
    #         A	     GCA	12	    41.38%           A   41.38
    #         A	     GCC	4	    13.79% ==>       C   66.67
    #         A	     GCT	12	    41.38%           E   64.44
    #         A	     GCG	1	     3.45%
    #         C	     TGC	4	    33.33%
    #         C 	     TGT	8	    66.67%
    #         E	     GAG	16	    35.56%
    #         E	     GAA	29	    64.44%

    ####### $$$$$$$$$$$$$$$$$$$$$$$$$$$$

    # compare all the condon's bias across all your samples using one way anova
    # /Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova/populationC_allsamples_bias.txt
    # df_allsamples_bias = pd.DataFrame(pd.read_csv("/Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova/populationC_allsamples_bias.txt", sep='\t'))
    #out_path = '/Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova/'


    ## here aim to combine all the populations
    df_allsamples_bias_all_t = pd.DataFrame()
    for file in os.listdir(out_path):

        if file.endswith('_allsamples_bias.txt'):
            file_path, file_name = os.path.split(file)
            prefix, middle, file_ext = file_name.split('_')

            df_allsamples_bias = pd.DataFrame(pd.read_csv(join(out_path, file), sep='\t'))
            df_allsamples_bias_t = df_allsamples_bias.T

            df_allsamples_bias_all_t = pd.concat([df_allsamples_bias_all_t,df_allsamples_bias_t],ignore_index=True,axis=1)

    #df_anova_all = pd.DataFrame(df_allsamples_bias_all_t.rename(columns=df.iloc[1]))
    df_anova_all = pd.DataFrame(df_allsamples_bias_all_t).to_csv(join(out_path,"allsamples_bias_t.txt"), sep="\t")

    #n = file_numbers(out_p)

 ###############################
    df_anova_all_noheader = df_anova_all.drop(df_anova_all.index[[0,1]]) # row is the sample, col is the bias
    df_anova_all_noheader[0]

    for i in range(0, 65):
        df_anova = df_anova_all_noheader[i]
        # get first each col
        df.rename(columns={x:y for x,y in zip(df.columns,range(0,len(df.columns)))})

        df_anova.columns = ['A','B','C','D']
        df_anova_1 = df_avova_1
        F, p = stats.f_oneway()


###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
























