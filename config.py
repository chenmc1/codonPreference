

import subprocess
import os
import pandas as pd


#def dictionary():

D1 = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
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


D2 = {}
for i in D1.keys():
    if D1[i] not in D2:
        D2[D1[i]] = {i: 0}
    else:
        D2[D1[i]][i] = 0


D3 = {i: 0 for i in D1.keys()}
D3['NNN'] = 0
#return D1, D2, D3

    # D2 = {
    # 'F': {'TTT': 0, 'TTC': 0}, 'S': {'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0},
    # 'Y': {'TAT': 0, 'TAC': 0}, 'C': {'TGT': 0, 'TGC': 0}, 'W': {'TGG': 0},
    # 'L': {'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'TTA': 0, 'TTG': 0},
    # 'P': {'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0},
    # 'I': {'ATT': 0, 'ATC': 0, 'ATA': 0}, 'M': {'ATG': 0},
    # 'R': {'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0, 'AGG': 0},
    # 'F': {'TTT': 0, 'TTC': 0}, 'F': {'TTT': 0, 'TTC': 0}, 'T': {'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0},
    # 'N': {'AAT': 0, 'AAC': 0}, 'K': {'AAA': 0, 'AAG': 0}, 'F': {'TTT': 0, 'TTC': 0}, 'Q': {'CAA': 0, 'CAG': 0},
    # 'V': {'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0}, 'A': {'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0},
    # 'G': {'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}, 'D': {'GAT': 0, 'GAC': 0}, 'E': {'GAA': 0, 'GAG': 0},
    # 'H': {'CAT': 0, 'CAC': 0},
    # 'STOP': {'TAA': 0, 'TAG': 0, 'TGA': 0, 'TGG': 0}}


    # D3 = {'TTT': 0, 'TTC': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'AGT': 0, 'AGC': 0,'TAT': 0, 'TAC': 0, 'TGT': 0,
    #       'TGC': 0, 'TGG': 0, 'CTT': 0, 'CTC': 0, 'CTA': 0,'CTG': 0, 'TTA': 0, 'TTG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0,
    #       'CCG': 0, 'ATT': 0,'ATC': 0, 'ATA': 0, 'ATG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 'AGA': 0,'AGG': 0,
    #       'TTT': 0, 'TTC': 0, 'TTT': 0, 'TTC': 0, 'ACT': 0, 'ACC': 0, 'ACA': 0,'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
    #       'TTT': 0, 'TTC': 0, 'CAA': 0, 'CAG': 0,'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GCT': 0, 'GCC': 0, 'GCA': 0,
    #       'GCG': 0,'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0,'CAT': 0, 'CAC': 0,
    #       'TAA': 0, 'TAG': 0, 'TGA': 0, 'TGG': 0, 'ACG': 0, "NNN": 0}





def check_dependencies(name):
    '''
    is_tool: aim to check if the required softwares that will been called on this main pipeline
    has been installed successfully.
    :param name: software name that can be call directly as command
    :return: False if the software is not installed or not in your path
    '''
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()

    except OSError as e:
        if e.errno == os.errno.ENOENT:

           return False
    return True


def checkHeader_readscountfile(each_readscountfile):
    '''
    check header for the reads count csv file, the headers for gene ids and reads count should be "Sequence_Name" and "Count"
    '''
    #each_readscountfile = "DGE0101_readcount_h.csv"
    Sequence_Name = 'Sequence_Name'
    Count = 'Count'

    input = pd.read_csv(each_readscountfile)
    headers = list(input.columns.values)

    if Sequence_Name in headers and Count in headers:

        #print('reads count file header passed check!')
        return True
    else:

        #print("File header Error: unable to find 'Sequence_Name' and 'Count' columns in your reads count file or you didn't name these two columns exactly as Sequence_Name and Count.")
        return False


def check_CDSfasta(each_cdsfile):


    # try:
    #     f = open(each_cdsfile)
    # except IOError:
    #     print('No such file: %s' % each_cdsfile)

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
            return False
        elif v[-3:].upper() not in ["TAA", "TAG", "TGA"]:
            print(
                "Warning: You have at least one transcript: %s that is not terminated by stop codon, be sure to estimate cds using TransDecoder" %k)
            return False
        elif Seq(v).translate().find("*") < len(v) / 3 - 1:
            print(
                "Warning: You have at least one transcript: %s that has a premature stop codon" %k)
            return False

        else:
            print("Input-cds-file quality pre-check passed!")

            return True


##










