__author__ = 'chenmingcui'


import os
import sys

filename = sys.argv[1]

def codonCount(filename):
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
       f = open(filename)
    except IOError:
       print('No such file: %s' %filename)

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
    f = open('count_result.txt','w')
    f.write(("AminoAcid\tCodon\tCount\n"))
    codon_count_pair = list(dict_2.values()) # values are in a list   [{'GAT': 78, 'GAC': 18}, {'GTG': 38, 'GTA': 20, 'GTT': 50, 'GTC': 24}]
    for n in range (0,len(codon_count_pair)):
           codons =list(codon_count_pair[n].keys())
           for m in range (0,len(codons)):
               count = codon_count_pair[n][codons[m]]
               f.write(dict_1[codons[m]] + '\t' + codons[m] + '\t' + str(count) + '\n')
    f.close()

    print('Counting finished, check your output file: count_reult.txt')

    return()


## now test it:

codonCount("sample.cds.fasta")


