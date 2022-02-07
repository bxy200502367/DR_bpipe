#!/bin/sh

MODULE_HOME="Module"
fastq_files=$1
NC=$2
PC=$3

cp -r $1 Sample

ls Sample > barcode_info.txt ##获得样本信息（barcode）

for i in $(cat barcode_info.txt)
do
{
gunzip -q Sample/$i/*
cat Sample/$i/* > $i.fastq
}&
done

wait

#############主要流程bpipe架构###########
bpipe run -r $MODULE_HOME/DR.bpipe *.fastq

#############删除冗余文件###########
rm -rf Sample
rm -f *.fastq

#############合并数据###########
for i in $(cat barcode_info.txt)
do
mkdir ${i}_analysis
mv ${i}_AMR_analysis ${i}_analysis
mv ${i}_SNP_analysis ${i}_analysis
mv ${i}_SGA_001_analysis ${i}_analysis
mv ${i}_SGA_002_analysis ${i}_analysis
mv ${i}_SGA_003_analysis ${i}_analysis
mv ${i}_MLST_analysis ${i}_analysis
done


#############整理结果###########
###SGA_001###
mkdir SGA_001_result_tmp
for i in $(cat barcode_info.txt)
do
cp ${i}_analysis/${i}_SGA_001_analysis/${i}.SGA_001.result.txt SGA_001_result_tmp/${i}.SGA_001.result.txt
done
cat SGA_001_result_tmp/* > SGA_001_result.txt
sed -i '1i Sample\tSGA_001_result\tDepth' SGA_001_result.txt
rm -rf SGA_001_result_tmp

###SGA_002###
mkdir SGA_002_result_tmp
for i in $(cat barcode_info.txt)
do
cp ${i}_analysis/${i}_SGA_002_analysis/${i}.SGA_002_result.txt SGA_002_result_tmp/${i}.SGA_002.result.txt
done
cat SGA_002_result_tmp/* > SGA_002.result.txt
awk 'NR%2 == 0 || NR == 1' SGA_002.result.txt > SGA_002_result.txt #提取第一行和偶数行
rm -f SGA_002.result.txt
rm -rf SGA_002_result_tmp

mkdir result
python $MODULE_HOME/sub/merge_result.py -i ${fastq_files}

#############删除冗余文件###########
rm -f SGA_001_result.txt
rm -f SGA_002_result.txt
rm -f barcode_info.txt

python $MODULE_HOME/sub/AMR_TF.py -p ${PC} -n ${NC}

rm -rf Module
