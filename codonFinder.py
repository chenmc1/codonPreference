__author__ = 'chenmingcui'

import os
from seq import read_fasta_file
import sys
import csv
#filename = sys.argv[1]
import os
os.chdir("/Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference")
# #def codon2aa(filename):
#     """
#     this module is to translate your codon to amino acid sequence from
#     your cds.fasta file.
#     workingDir:/Users/chenmingcui/Documents/A-PhD_research_projects/PhD_Project/codonPreference
#     """
## generate dic_1 : 'codon 2 amino acids'
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

## generate dic_2 : 'amino acids 2 codons count'
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

# check if the dict_1 and dict_2 are matched and completed
aa_dict1 = list(dict_1.values())
aa_dict2 = list(dict_2.keys())

# check file existance
try:
   f = open("sample.cds.fasta")
except IOError:
    print('No such file!')
    #print('File %s does not exit!!' % filename)
# read a fasta file
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

# access to the transcript and count the codons

seq =list(seqs.values())   # convert dictionary values into list otherwise cannot be indexed
for i in range (0,len(seq)):
    each_seq = seq[i]
    for s in range (0,len(each_seq),3):
    #codon_count_dic = dict_2[dict_1[each_seq[s:s+3].upper()]]
        codon = each_seq[s:s+3].upper()
        aa = dict_1[codon]
        codon_count_dict = dict_2[aa]
        codon_count_dict[codon] += 1

for aa, codon_count in dict_2.items():
    print (aa,codon_count)


#result_list =[]
#aa_result =list(dict_2.keys())            # keys are in a list
#print("the genome of the file : %s" + %filename)
f = open('count_result.txt','w')
print("|AminoAcid|Codon|Count|")
f.write(("AminoAcid\tCodon\tCount\n"))
codon_count_pair = list(dict_2.values()) # values are in a list   [{'GAT': 78, 'GAC': 18}, {'GTG': 38, 'GTA': 20, 'GTT': 50, 'GTC': 24}]
for n in range (0,len(codon_count_pair)):
       codons =list(codon_count_pair[n].keys())
       for m in range (0,len(codons)):
           count = codon_count_pair[n][codons[m]]
           #each_row = dict_1[codons[m]] + ',' + codons[m] + ',' + str(count)
           #each_row = each_row.split(',')
           #print("|" + dict_1[codons[m]] + " "*(9-len(dict_1[codons[m]])) + "|" + codons[m] + " "*(5-len(codons[m])) + "|" + str(count) + " "*(5-len(str(count))) + "|")
#            result = print("|" + dict_1[codons[m]] + " "*(9-len(dict_1[codons[m]])) + "|" + codons[m] + " "*(5-len(codons[m])) + "|" + str(count) + " "*(5-len(str(count))) + "|")
#
#            with open('outfile.txt','w') as outfile:
#               outfile.write(result)
#print("|" + dict_1[codons[m]] + " "*(9-len(dict_1[codons[m]])) + "|" + codons[m] + " "*(5-len(codons[m])) + "|" + str(count) + " "*(5-len(str(count))) + "|")
# outfile.close()
           print( dict_1[codons[m]] + codons[m]  + str(count))
           f.write(dict_1[codons[m]] + '\t' + codons[m] + '\t' + str(count) + '\n')
f.close()


           f = open('count_result.txt','w')
           result = "|" + dict_1[codons[m]] + " "*(9-len(dict_1[codons[m]])) + "|" + codons[m] + " "*(5-len(codons[m])) + "|" + str(count) + " "*(5-len(str(count))) + "|"
           f.write(result)
           f.close()



          # each_row = dict_1[codons[m]] + "\t" + codons[m] + "\t" + str.count + "\n"
          # #each_row = each_row.split()
          # table = table + each_row
          #
          # with open('result.csv','w') as csvfile:
          #      fieldnames = ['Amino Acid', 'Codon', 'Count']
          #      writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
          #      writer.writeheader()
          #      writer.writerow(dict_1[codons[m]] + "," + codons[m] + str(count) + "\n")
          #
          #
