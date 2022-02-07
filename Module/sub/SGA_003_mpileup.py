"""
@author = 'YuanXu'
@date = '2022/01/26'
功能：三代靶向耐药流程--SGA_003模块
description：kpc基因研究衍生，先从mpileup文件拟合生成一个fasta(包括点突变、Indel以及去除前后未测序到的片段)
"""

import argparse
import os
import pandas as pd
from re import match

def myParser():
    parser = argparse.ArgumentParser(description = 'SNP module')
    parser.add_argument('-i', '--input_fastq', help = "输入测序原始数据，fastq格式文件")
    parser.add_argument('-d', '--Indel_info', help = "输入插入缺失的信息，vcf格式")
    args = parser.parse_args()
    return args

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
        for i in Gene_name_list:
            fit_sequence = ''.join(grouped2.get_group(i)['Alt'].tolist())
            output_mpileup_file_name = output_path + '/' + abstract_barcode_name(input_fastq_file) + "." + i + ".fit.fasta"
            with open(output_mpileup_file_name,'w') as out_fasta_fit:
                out_fasta_fit.write('>' + abstract_barcode_name(input_fastq_file) + '_' + i + '\n')
                out_fasta_fit.write(fit_sequence + '\n')

def abstract_barcode_name(input_fastq_file):
    dirname, basename = os.path.split(input_fastq_file)
    barcode_name = basename.split('.')[0]
    return barcode_name





if __name__ == '__main__':
    args = myParser()
    
    input_fastq_file = args.input_fastq
    barcode_name = abstract_barcode_name(input_fastq_file)
    mpileup_file = f"{barcode_name}.SGA_003.mpileup.txt"
    
    with open(mpileup_file, 'r') as f:
        mpileup_hits_list = f.read().rstrip().split("\n") #列表，一行为一元素
    mpileup_hits = MpileupHits(mpileup_hits_list) #把mpileup结果解析成MpileuppHits对象
    mpileup_hits.write_mpileup(input_fastq_file, 'Mpileup_count')
    
    Indel_info_vcf = args.Indel_info
    with open(Indel_info_vcf, 'r') as f:
        Indel_info = f.read().rstrip().split("\n")
    del Indel_info[0]
    
    Indel_info_dict = {}
    for i in Indel_info:
        chrom,position,Ref,Var = tuple(i.split("\t"))[0:4]
        Indel_info_dict[(chrom,position)] = Var
    
    df = mpileup_hits.pd_arrange.copy(deep = True)
    for k,v in Indel_info_dict.items(): #Indel替换
        gene_name = k[0]
        position = k[1]
        Indel = v
        if Indel.startswith('-'):
            Indel_length = len(Indel) - 1
            for i in range(Indel_length):
                del_position = int(position) + i + 1
                row_index = df[(df['Gene_name'] == k[0]) & (df['Pos'] == str(del_position))].index[0]
                df.loc[row_index, 'Alt'] = '-'
                df.loc[row_index, 'Genotype'] = '-'
        if Indel.startswith('+'):
            insert_position = int(position) + 1
            row_index = df[(df['Gene_name'] == k[0]) & (df['Pos'] == str(insert_position))].index[0]
            df.loc[row_index, 'Alt'] = Indel
            df.loc[row_index, 'Genotype'] = Indel
            
    df2 = df.copy(deep = True)
    Gene_name_list = list(set(df2['Gene_name']))
    grouped2 = df2.groupby('Gene_name')
    for i in Gene_name_list:
        fit_sequence = ''.join(grouped2.get_group(i)['Alt'].tolist())
        fit_sequence = fit_sequence.replace('-','')
        fit_sequence = fit_sequence.replace('+','')
        fit_sequence.strip('N')
        output_fit_fasta = 'Fit_fasta/' + abstract_barcode_name(input_fastq_file) + "." + i + ".fit.fasta"
        with open(output_fit_fasta,'w') as f:
            f.write('>' + abstract_barcode_name(input_fastq_file) + '_' + i + '\n')
            f.write(fit_sequence + '\n')
