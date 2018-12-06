import os
from os.path import join
import pandas as pd


def combine_samples(in_path, out_path):
    """
    since you have multiple samples for each population or treatment, 'combine_sample' is to combine the bias.txt files of each population/treatment into one file.

              population_A/TrtA         population_B/TrtB         populationC/TrtC
              sample_A1                 sample_B1                 sample_C1
              sample_A2                 sample_B2                 sample_C2
              sample_A3                 sample_B3                 sample_C3
                    :
                    :
    For example, here is to comebine your sample_A1,A2,A3 into population_A... comebine your sample_B1,B2,B3 into population_B...

    Input:
           Note that here we just  process pop-A_sample-1_countBias_result.txt files that from same population/trts, so if you have 3 populations you need to call
           this script 3 times. you can put it in a loop to run for multiple populations.

           - input path for these files
           - output path

    :return: populationA_condonbias.txt
    """
    #in_path='/Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova'
    #out_path='/Users/chenmingcui/Documents/PhD_work/trivial_scripts/test_anova'
    # break path and fileanme into prefix
    all_file_names = []
    all_pop_names = []
    #[file for file in os.listdir(in_dir) if file.endswith('bias_count_result.txt')]
    for file in os.listdir(in_path):
        if file.endswith('bias_count_result.txt'):
            file_path, file_name  = os.path.split(file)
            prefix, middle, file_ext = file_name.split('.')
            population_name, sample_name = prefix.split('-')
            all_file_names.append(file_name)
            all_pop_names.append(population_name)

    all_pop_names = sorted(set(all_pop_names))

    dict_all = dict([(key, []) for key in all_pop_names])

    # summary the input files into a dictionary
    for i in range(0,len(all_file_names)):
        for key in dict_all:
            if all_file_names[i][0:11] == key:
                dict_all[key].append(all_file_names[i])

    # update dictionary like below:

    # {'populationA': ['populationA-sampleA1.bias_count_result.txt',
    #                  'populationA-sampleA2.bias_count_result.txt',
    #                  'populationA-sampleA3.bias_count_result.txt'],
    #  'populationB': ['populationB-sampleB1.bias_count_result.txt',
    #                  'populationB-sampleB2.bias_count_result.txt',
    #                  'populationB-sampleB3.bias_count_result.txt'],
    #  'populationC': ['populationC-sampleC1.bias_count_result.txt',
    #                  'populationC-sampleC2.bias_count_result.txt',
    #                  'populationC-sampleC3.bias_count_result.txt']}

    for key in dict_all:
        each_file_list = dict_all.get(key)
        #df_codonbias = pd.DataFrame()
        #print(each_file_list)
        appended_data = []
        for each_file in each_file_list:
            data = pd.read_csv(join(in_path,each_file),sep='\t')
            appended_data.append(data)
        appended_data = pd.concat(appended_data, ignore_index=True, axis=1) # combine all files in a list into one df

        print("with "+key+"\n",appended_data)

        appended_data.to_csv(join(out_path,key+'_combined_codonbias.txt'), sep='\t')

        print(key+" write into file")















