"""
@author = 'YuanXu'
@date = '2022/01/20'
功能：三代靶向耐药流程--SNP模模块
description：通过mpileup文件查找点突变结果并转化成氨基酸
"""

import argparse
import re
from re import split
from re import match
import os
import numpy as np
import pandas as pd
import subprocess as sp

from Bio import SeqIO

def myParser():
    parser = argparse.ArgumentParser(description = 'SNP module')
    parser.add_argument('-i', '--input_fastq', help = "输入测序原始数据，fastq格式文件或者文件夹")
    args = parser.parse_args()
    return args

class ReferenceFasta(object):
    def __init__(self,reference_fasta_file):
    
        with open(reference_fasta_file,"r") as ref_fa:
            self.ref_fa_dict = {}
            for line in ref_fa:
                line= line.strip() #strip()方法移除字符串头尾指定的字符（默认为空格或换行符）,不会删除中间部分字符
                if line == '':
                    continue
                if line.startswith('>'):
                    read_name = line.lstrip('>')
                    read_name_dict = re.sub('\..*', '', read_name)
                    '''把>标签行整理成bed文件中read标签的格式(比如fasta文件中>标签行：'gi|41406098|ref|NC_002944.2| Mycobacterium avium subsp. paratuberculosis K-10, complete sequence'（ncbi标准格式）；bed文件中第一列read标签：如NC_002944.2（只需其中一部分）)，运用自己构建的正则表达提取。例句只提取第一个空格前的'''
                    self.ref_fa_dict[read_name] = ''
                else:
                    self.ref_fa_dict[read_name] += line
                    
        self.ref_fa_names = pd.DataFrame(self.ref_fa_dict.keys(), columns = ["Gene_names"], dtype = object)
        
        ref_fa_names = self.ref_fa_names.copy(deep = True)
        ref_fa_names.loc[ref_fa_names.shape[0]] = 'Total'
        ref_fa_names.loc[ref_fa_names.shape[0]] = 'Total-phage'
        ref_fa_names['Depth'] = int(0)
        ref_fa_names['Ave_Length'] = ''
        ref_fa_names['Ave_Score'] = ''
        ref_fa_names['Ave_Identity'] = ''
        self.ref_fa_count_initialization = ref_fa_names
        
        ref_fa_names_SNP = self.ref_fa_names.copy(deep = True)
        ref_fa_names_SNP['Depth'] = int(0)
        self.ref_fa_count_initialization_SNP = ref_fa_names_SNP
        
        self.protein_dict = {}
        for k,v in self.ref_fa_dict.items():
            self.protein_dict[k] = translate_DNA(v)
            
    def count_result_depth_SNP(self,filter_result_count_pd):
        result_depth_count_SNP = self.ref_fa_count_initialization_SNP.copy(deep=True)
        for i in filter_result_count_pd['Gene_names'].to_list():
            result_depth_count_SNP.loc[result_depth_count_SNP['Gene_names'] == i, 'Depth'] = filter_result_count_pd.loc[i, 'Depth']
        return result_depth_count_SNP


def extractCigarSeq(sequence, ref, phred, mapq):
    letters = dict([('A', list()), ('C', list()), ('G', list()), ('T', list()), ('N', list()), ('Insert', list()), ('Del', list())])
    lettersKeys = letters.keys()
    sequence = sequence.upper() #转换成大写，比对到正向链上大写，比对到互补链上用小写
    
    asciiOffset = 33 #质量都是Phred3的值
    positionSeq = 0
    positionPhred = 0
    positionMapq = 0
    
    while(positionSeq < len(sequence)):
        currentData = sequence[positionSeq]
        if currentData == "," or currentData == ".": #比对匹配的情况
            phredVal = ord(phred[positionPhred]) - asciiOffset
            mapqVal  = ord(mapq[positionMapq]) - asciiOffset
            letters[ref].append((phredVal, mapqVal))
            positionSeq += 1
            positionPhred += 1
            positionMapq += 1
            
        elif currentData in lettersKeys: #单位点不匹配的情况
            phredVal = ord(phred[positionPhred]) - asciiOffset
            mapqVal  = ord(mapq[positionMapq]) - asciiOffset
            letters[currentData].append((phredVal, mapqVal))
            positionSeq += 1
            positionPhred += 1
            positionMapq += 1
            
        elif currentData == "+": #存在插入缺失的情况
            letters['Insert'].append((-1, -1))
            res = match(r"[+](\d+)", sequence[positionSeq:len(sequence)])
            indelLength = res.groups()[0]
            positionSeq = positionSeq + 1 + len(indelLength) + int(indelLength)
            
        elif currentData == "-": #存在插入缺失的情况
            letters['Del'].append((-1, -1))
            res = match(r"[-](\d+)", sequence[positionSeq:len(sequence)])
            indelLength = res.groups()[0]
            positionSeq = positionSeq + 1 + len(indelLength) + int(indelLength)
            
        elif currentData == ">" or currentData == "<":
            positionSeq += 1
            positionPhred += 1
            positionMapq += 1
            
        elif currentData == "$": #$代表一个read的结束，该符号修饰的是其前面的碱基。
            positionSeq += 1
            
        elif currentData == "*": #"*"代表模糊碱基
            positionSeq += 1
            
        elif currentData == "^": # "^"代表匹配的碱基是一个read的开始；’^'后面紧跟的ascii码减去33代表比对质量；这两个符号修饰的是后面的碱基，其后紧跟的碱基(.,ATCGatcgNn)代表该read的第一个碱基。
            positionSeq += 2
        else:
            print("Problem with letter: " + currentData)
            print("in this cigar string: " + sequence)
            sys.exit(2)
            
    return(letters)

class MpileupHit(object):
    def __init__(self,line):
        self.line = line.rstrip()
        fields = line.rstrip().split("\t")
        self.gene_name = fields[0]
        self.pos = int(fields[1])
        self.ref = fields[2]
        self.depth = int(fields[3])
        self.map_info = fields[4]
        self.phred = fields[5]
        self.mapq = fields[6]
        
        self.letters = extractCigarSeq(self.map_info, self.ref, self.phred, self.mapq)
        self.lettersA = len(self.letters['A'])
        self.lettersC = len(self.letters['C'])
        self.lettersG = len(self.letters['G'])
        self.lettersT = len(self.letters['T'])
        self.lettersN = len(self.letters['N'])
        self.lettersInsert = len(self.letters['Insert'])
        self.lettersDel = len(self.letters['Del'])

class MpileupHits(object):
    def __init__(self, mpileup_hits_list):
        self.mpileup_hits = []
        self.mpileup_hits_pd_list = []
        self.mpileup_lettersA = []
        self.mpileup_lettersC = []
        self.mpileup_lettersG = []
        self.mpileup_lettersT = []
        self.mpileup_lettersN = []
        self.mpileup_lettersInsert = []
        self.mpileup_lettersDel = []
        for i in mpileup_hits_list:
            self.mpileup_hits.append(MpileupHit(i))
            self.mpileup_lettersA.append(MpileupHit(i).lettersA)
            self.mpileup_lettersC.append(MpileupHit(i).lettersC)
            self.mpileup_lettersG.append(MpileupHit(i).lettersG)
            self.mpileup_lettersT.append(MpileupHit(i).lettersT)
            self.mpileup_lettersN.append(MpileupHit(i).lettersN)
            self.mpileup_lettersInsert.append(MpileupHit(i).lettersInsert)
            self.mpileup_lettersDel.append(MpileupHit(i).lettersDel)
            self.mpileup_hits_pd_list.append(i.split("\t"))
        
        self.mpileup_hits_pd = pd.DataFrame(self.mpileup_hits_pd_list)
        self.mpileup_hits_pd['A'] = self.mpileup_lettersA
        self.mpileup_hits_pd['C'] = self.mpileup_lettersC
        self.mpileup_hits_pd['G'] = self.mpileup_lettersG
        self.mpileup_hits_pd['T'] = self.mpileup_lettersT
        self.mpileup_hits_pd['N'] = self.mpileup_lettersN
        self.mpileup_hits_pd['Insert'] = self.mpileup_lettersInsert
        self.mpileup_hits_pd['Del'] = self.mpileup_lettersDel
        
        columns_names = ['Gene_name', 'Pos', 'Ref', 'Depth', 'map_info', 'phred', 'mapq', 'A', 'C', 'G', 'T', 'N', 'Insert', 'Del']
        self.mpileup_hits_pd.columns = columns_names
        del self.mpileup_hits_pd['map_info']
        del self.mpileup_hits_pd['phred']
        del self.mpileup_hits_pd['mapq']
        
        self.pd_arrange = self.mpileup_hits_pd.copy(deep = True)
        data_types_dict = {'Depth': int, 'A': int, 'C': int, 'G': int, 'G': int, 'N': int, 'Insert': int, 'Del': int} #转换数据格式
        self.pd_arrange['Depth'] = self.pd_arrange['A'] + self.pd_arrange['C'] + self.pd_arrange['G'] + self.pd_arrange['T']
        self.pd_arrange = self.pd_arrange.astype(data_types_dict)
        self.pd_arrange['Ref_num'] = [self.pd_arrange.loc[i,self.pd_arrange.loc[i,'Ref']] for i in range(self.pd_arrange.shape[0])]
        self.pd_arrange['Ref_fre'] = round(self.pd_arrange['Ref_num'] / self.pd_arrange['Depth']*100,2).astype(str)+'%'
        self.pd_arrange['Alt'] = self.pd_arrange[['A','C','G','T']].idxmax(1) #ACGT中数值最大的索引
        del self.pd_arrange['Ref_num']
        self.pd_arrange['Genotype'] = self.pd_arrange['Ref'] + '>' + self.pd_arrange['Alt']
        self.pd_arrange.loc[self.pd_arrange['Depth'] == 0, 'Alt'] = 'N'
        self.pd_arrange.loc[self.pd_arrange['Depth'] == 0, 'Genotype'] = 'N'
        
    def write_mpileup(self, input_fastq_file, output_path):
        df = self.pd_arrange.copy(deep = True)
        Gene_name_list = list(set(df['Gene_name']))
        grouped = df.groupby('Gene_name')
        for i in Gene_name_list:
            output_mpileup_file_name = output_path + '/' + abstract_barcode_name(input_fastq_file) + "." + i + ".mpileup.xlsx"
            grouped.get_group(i).to_excel(output_mpileup_file_name, index=False)
            
    def write_fasta_fit(self, input_fastq_file, output_path):
        df2 = self.pd_arrange.copy(deep = True)
        Gene_name_list = list(set(df2['Gene_name']))
        grouped2 = df2.groupby('Gene_name')
        protein_dict = {}
        
        for i in Gene_name_list:
            fit_sequence = ''.join(grouped2.get_group(i)['Alt'].tolist())
            
            protein = translate_DNA(fit_sequence)
            protein_dict[i] = protein
            
            output_mpileup_file_name = output_path + '/' + abstract_barcode_name(input_fastq_file) + "." + i + ".fit.fasta"
            with open(output_mpileup_file_name,'w') as out_fasta_fit:
                out_fasta_fit.write('>' + abstract_barcode_name(input_fastq_file) + '_' + i + '\n')
                out_fasta_fit.write(fit_sequence + '\n')
            out_fasta_fit.close()
        return protein_dict

def abstract_barcode_name(input_fastq_file):
    dirname, basename = os.path.split(input_fastq_file)
    barcode_name = basename.split('.')[0]
    return barcode_name

def translate_DNA(dna):
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    
    amino_acid_sequence = ""
    for start in range(0,len(dna) - 2, 3):
        stop = start + 3
        codon = dna[start:stop]
        aa = gencode.get(codon.upper(),'X') #当指定键的值不存在时，返回X
        amino_acid_sequence = amino_acid_sequence + aa
    return(amino_acid_sequence)

def protein_mut(protein1, protein2):
    mut = []
    for i in range(len(protein1)):
        if protein2[i] != "X":
            if protein1[i] != protein2[i]:
                mut.append(protein1[i] + str(i+1) + protein2[i])
    return(mut)

def protein_mut_dict(protein_dict1, protein_dict2):
    mut_dict = {}
    for gene_name,protein in protein_dict2.items():
        if gene_name in protein_dict1.keys():
            mut_dict[gene_name] = protein_mut(protein_dict1[gene_name], protein)
    return mut_dict

def SNP_result_pd(protein_mut_info, ref_fa):
    count_df = ref_fa.copy(deep = True)
    count_df.index = count_df['Gene_names']
    count_df['mut_info'] = ''
    for gene_name,protein in protein_mut_info.items():
        count_df.loc[gene_name,'mut_info'] = ','.join(protein_mut_info[gene_name])
    del count_df['Gene_names']
    del count_df['Depth']
    return count_df

def run_SNP_module(input_fastq_file):
    barcode_name = abstract_barcode_name(input_fastq_file) #得到barcode/sample name
    reference_file = f"../Module/Reference/TMP/{barcode_name}_analysis/SNP/SNP.fasta"
    ref_fa = ReferenceFasta(reference_file) #ReferenceFasta对象，SNP模块需要ref_fa_names属性
    
    mpileup_file = f"{barcode_name}.SNP.mpileup.txt"
    with open(mpileup_file, 'r') as f:
        mpileup_hits_list = f.read().rstrip().split("\n") #列表，一行为一元素
    mpileup_hits = MpileupHits(mpileup_hits_list) #把mpileup结果解析成MpileupHits对象
    mpileup_hits.write_mpileup(input_fastq_file, 'Mpileup_count')
    
    protein_dict2 = mpileup_hits.write_fasta_fit(input_fastq_file, 'Fit_fasta') #返回拟合的核苷酸序列
    protein_dict1 = ref_fa.protein_dict
    protein_mut_info = protein_mut_dict(protein_dict1, protein_dict2) #氨基酸突变列表
    ref_fa_count_initialization_SNP = ref_fa.ref_fa_count_initialization_SNP
    SNP_result = SNP_result_pd(protein_mut_info, ref_fa_count_initialization_SNP) #参数protein_mut_info（字典）；ref_fa（ReferenceFasta对象的pd）,结果是pd
    SNP_result.rename(columns = {'mut_info' : barcode_name}, inplace=True) #把列名Depth改成barcode名字
    output_result_file = barcode_name + '_SNP_result.txt'
    SNP_result.to_csv(output_result_file, index=True, sep='\t')
    return SNP_result

if __name__ == '__main__':
    args = myParser()
    input_fastq_file = args.input_fastq
    
    SNP_result = run_SNP_module(input_fastq_file)
