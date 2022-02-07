sample_name=`echo $1 | cut -d. -f1`
MODULE_HOME=$2
TMP_reference=${MODULE_HOME}/Reference/TMP/${sample_name}_analysis

mkdir ${sample_name}_SNP_analysis
cd ${sample_name}_SNP_analysis
mkdir -p ${TMP_reference}/SNP

###################要分析点突变的基因##############################
cat ${MODULE_HOME}/Reference/SNP/* | dos2unix | seqkit seq -u > ${TMP_reference}/SNP/SNP.fasta #reference文件夹下的所有参考序列合并
java -jar /home/chensy/software/picard.jar CreateSequenceDictionary -REFERENCE ${TMP_reference}/SNP/SNP.fasta -OUTPUT ${TMP_reference}/SNP/SNP.dict #使用picard软件生成dict
lastdb ${TMP_reference}/SNP/SNP.db ${TMP_reference}/SNP/SNP.fasta #根据参考基因构建数据库

grep ">" ${TMP_reference}/SNP/SNP.fasta | sed "s/>//g" > ${TMP_reference}/SNP/reference_info_SNP.txt
echo "Gene_names" > ${TMP_reference}/SNP/reads_stat.txt
cat ${TMP_reference}/SNP/reads_stat.txt ${TMP_reference}/SNP/reference_info_SNP.txt > ${TMP_reference}/SNP/reads_stat_info.txt
##############################################################################

lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 ${TMP_reference}/SNP/SNP.db ../$1 | maf-convert -f ${TMP_reference}/SNP/SNP.dict sam -r 'ID:[id] PL:[nanopore] SM:[sample]' - > ${sample_name}.SNP.sam #last比对并且maf格式文件转换成sam格式文件
lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 -f BlastTab+ ${TMP_reference}/SNP/SNP.db ../$1 | grep -v ^# - > ${sample_name}.NanoFilt.BlastTab+ #比对生成blasttab+格式文件

grep "^@" ${sample_name}.SNP.sam > ${sample_name}.tmp1.sam #提取表头
grep -v "^@" ${sample_name}.SNP.sam > ${sample_name}.tmp2.sam #去除表头

if test -s ${sample_name}.tmp2.sam; then
Rscript ${MODULE_HOME}/sub/align_filter.R ${sample_name}.NanoFilt.BlastTab+ ${sample_name}.tmp2.sam ${sample_name}.NanoFilt.BlastTab+.txt ${sample_name}.tmp2.sam_filtered
else
touch ${sample_name}.NanoFilt.BlastTab+.txt
touch ${sample_name}.tmp2.sam_filtered
fi

cat ${sample_name}.tmp1.sam ${sample_name}.tmp2.sam_filtered > ${sample_name}.SNP.sam_filtered #生成过滤完的sam文件

rm -f ${sample_name}.tmp1.sam
rm -f ${sample_name}.tmp2.sam
rm -f ${sample_name}.tmp2.sam_filtered
rm -f ${sample_name}.SNP.sam
rm -f ${sample_name}.NanoFilt.BlastTab+

samtools view -bS ${sample_name}.SNP.sam_filtered | samtools sort - > ${sample_name}.SNP.sorted.bam #sam文件转化成bam文件并排序
samtools index ${sample_name}.SNP.sorted.bam #bam文件索引
samtools mpileup -a -s -f ${TMP_reference}/SNP/SNP.fasta ${sample_name}.SNP.sorted.bam > ${sample_name}.SNP.mpileup.txt

mkdir Fit_fasta
mkdir Mpileup_count
if test -s ${sample_name}.SNP.mpileup.txt; then
python ${MODULE_HOME}/sub/SNP_mpileup.py -i ../$1
else
touch tmp_SNP_result.txt
echo "$sample_name" > tmp_SNP_result.txt
paste ${TMP_reference}/SNP/reads_stat_info.txt tmp_SNP_result.txt > ${sample_name}_SNP_result.txt
rm -f tmp_SNP_result.txt
fi

rm -f ${sample_name}.SNP.sam_filtered
rm -f ${sample_name}.SNP.sorted.bam
rm -f ${sample_name}.SNP.sorted.bam.bai
rm -f ${sample_name}.SNP.mpileup.txt

###统计测序深度

mkdir abstract_fastq
mkdir abstract_fasta

touch ${TMP_reference}/SNP/${sample_name}.stat.txt
echo "Depth" >> ${TMP_reference}/SNP/${sample_name}.stat.txt

for i in $(cat ${TMP_reference}/SNP/reference_info_SNP.txt) #统计样本比对到的reads数
do
grep -w ${i} ${sample_name}.NanoFilt.BlastTab+.txt | awk '{print $1}' > ${sample_name}.${i}.abstract.txt #提取比对上的序列名

if test -s ${sample_name}.${i}.abstract.txt; then
seqkit grep -f ${sample_name}.${i}.abstract.txt ../$1 > ${sample_name}.${i}.abstract.fastq #提取比对上的序列
seqkit fq2fa ${sample_name}.${i}.abstract.fastq -o ${sample_name}.${i}.abstract.fasta #fastq转化成fasta
mv ${sample_name}.${i}.abstract.fastq abstract_fastq
mv ${sample_name}.${i}.abstract.fasta abstract_fasta
fi

let count=`cat ${sample_name}.${i}.abstract.txt | wc -l` #得到比对到的reads数
echo "$count" >> ${TMP_reference}/SNP/${sample_name}.stat.txt

rm -f ${sample_name}.${i}.abstract.txt
done

paste ${sample_name}_SNP_result.txt ${TMP_reference}/SNP/${sample_name}.stat.txt > ${sample_name}.SNP.result.txt

rm -f ${sample_name}.NanoFilt.BlastTab+.txt
rm -f ${sample_name}_SNP_result.txt
