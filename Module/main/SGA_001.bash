sample_name=`echo $1 | cut -d. -f1`
MODULE_HOME=$2
TMP_reference=${MODULE_HOME}/Reference/TMP/${sample_name}_analysis

mkdir ${sample_name}_SGA_001_analysis
cd ${sample_name}_SGA_001_analysis
mkdir -p ${TMP_reference}/SGA_001 #创建存放reference的临时文件

###################要分析点突变的基因##############################
cat ${MODULE_HOME}/Reference/SGA_001/* | dos2unix | seqkit seq -u > ${TMP_reference}/SGA_001/SGA_001.fasta #SGA_001下所有参考序列合并
java -jar /home/chensy/software/picard.jar CreateSequenceDictionary -REFERENCE ${TMP_reference}/SGA_001/SGA_001.fasta -OUTPUT ${TMP_reference}/SGA_001/SGA_001.dict #使用picard软件生成dict
lastdb ${TMP_reference}/SGA_001/SGA_001.db ${TMP_reference}/SGA_001/SGA_001.fasta #根据参考基因构建数据库
grep ">" ${TMP_reference}/SGA_001/SGA_001.fasta | sed "s/>//g" > ${TMP_reference}/SGA_001/reference_info_SGA_001.txt
echo "Gene_names" > ${TMP_reference}/SGA_001/reads_stat.txt
cat ${TMP_reference}/SGA_001/reads_stat.txt ${TMP_reference}/SGA_001/reference_info_SGA_001.txt > ${TMP_reference}/SGA_001/reads_stat_info.txt
##############################################################################

lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 ${TMP_reference}/SGA_001/SGA_001.db ../$1 | maf-convert -f ${TMP_reference}/SGA_001/SGA_001.dict sam -r 'ID:[id] PL:[nanopore] SM:[sample]' - > ${sample_name}.SGA_001.sam #last比对并且maf格式文件转换成sam格式文件
lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 -f BlastTab+ ${TMP_reference}/SGA_001/SGA_001.db ../$1 | grep -v ^# - > ${sample_name}.NanoFilt.BlastTab+ #比对生成blasttab+格式文件

grep "^@" ${sample_name}.SGA_001.sam > ${sample_name}.tmp1.sam #提取表头
grep -v "^@" ${sample_name}.SGA_001.sam > ${sample_name}.tmp2.sam #去除表头

if test -s ${sample_name}.tmp2.sam; then
Rscript ${MODULE_HOME}/sub/align_filter.R ${sample_name}.NanoFilt.BlastTab+ ${sample_name}.tmp2.sam ${sample_name}.NanoFilt.BlastTab+.txt ${sample_name}.tmp2.sam_filtered
else
touch ${sample_name}.NanoFilt.BlastTab+.txt
touch ${sample_name}.tmp2.sam_filtered
fi

cat ${sample_name}.tmp1.sam ${sample_name}.tmp2.sam_filtered > ${sample_name}.SGA_001.sam_filtered

rm -f ${sample_name}.tmp1.sam
rm -f ${sample_name}.tmp2.sam
rm -f ${sample_name}.tmp2.sam_filtered
rm -f ${sample_name}.SGA_001.sam
rm -f ${sample_name}.NanoFilt.BlastTab+
#rm -f ${sample_name}.NanoFilt.BlastTab+.txt

samtools view -bS ${sample_name}.SGA_001.sam_filtered | samtools sort - > ${sample_name}.SGA_001.sorted.bam #sam文件转化成bam文件并排序
samtools index ${sample_name}.SGA_001.sorted.bam #bam文件索引
samtools mpileup -a -s -f ${TMP_reference}/SGA_001/SGA_001.fasta ${sample_name}.SGA_001.sorted.bam > ${sample_name}.SGA_001.mpileup.txt

mkdir Fit_fasta
mkdir Mpileup_count

if test -s ${sample_name}.SGA_001.mpileup.txt; then
samtools mpileup -f ${TMP_reference}/SGA_001/SGA_001.fasta ${sample_name}.SGA_001.sorted.bam | java -jar ${MODULE_HOME}/sub/VarScan.v2.3.9.jar mpileup2indel --min-var-freq 0.5 –output-vcf 1 > ${sample_name}.SGA_001.InDel.vcf #Varscan生成插入缺失
sed -i '1d' ${sample_name}.SGA_001.InDel.vcf #删除第一行
python ${MODULE_HOME}/sub/SGA_001_mpileup.py -i ../$1
else
touch ${sample_name}.SGA_001.InDel.vcf
touch ${sample_name}_SGA_001_result.txt
fi

rm -f ${sample_name}.SGA_001.sam_filtered
rm -f ${sample_name}.SGA_001.sorted.bam
rm -f ${sample_name}.SGA_001.sorted.bam.bai
rm -f ${sample_name}.SGA_001.mpileup.txt

#################################统计测序深度

mkdir abstract_fastq
mkdir abstract_fasta

touch ${TMP_reference}/SGA_001/${sample_name}.stat.txt

echo "sample" > ${TMP_reference}/SGA_001/reads_stat.txt
cat ${TMP_reference}/SGA_001/reads_stat.txt ${TMP_reference}/SGA_001/reference_info_SGA_001.txt > ${TMP_reference}/SGA_001/reads_stat_info.txt

for i in $(cat ${TMP_reference}/SGA_001/reference_info_SGA_001.txt) #统计样本比对到的reads数
do
grep -w ${i} ${sample_name}.NanoFilt.BlastTab+.txt | awk '{print $1}' > ${sample_name}.${i}.abstract.txt #提取比对上的序列名

if test -s ${sample_name}.${i}.abstract.txt; then
seqkit grep -f ${sample_name}.${i}.abstract.txt ../$1 > ${sample_name}.${i}.abstract.fastq #提取比对上的序列
seqkit fq2fa ${sample_name}.${i}.abstract.fastq -o ${sample_name}.${i}.abstract.fasta #fastq转化成fasta
##合并结果
mv ${sample_name}.${i}.abstract.fastq abstract_fastq
mv ${sample_name}.${i}.abstract.fasta

fi

let count=`cat ${sample_name}.${i}.abstract.txt | wc -l` #得到比对到的reads数
echo "$count" >> ${TMP_reference}/SGA_001/${sample_name}.stat.txt
touch PA-4_result.txt

if [ $count -le 10 ] ; then
echo 'LD' > PA-4_result.txt
else
############################PA-4结果判断
if test -s ${sample_name}.SGA_001.InDel.vcf; then
echo "SM" > PA-4_result.txt
else
       if test -s ${sample_name}_SGA_001_result.txt; then
       echo "WM" > PA-4_result.txt
       else
	   echo "W" > PA-4_result.txt
	   fi
fi
fi
rm -f ${sample_name}.${i}.abstract.txt
done
paste ${TMP_reference}/SGA_001/reads_stat_info.txt ${TMP_reference}/SGA_001/${sample_name}.stat.txt > ${sample_name}.SGA_001.depth_result.txt
echo ${sample_name} >> tmp.sample_name.txt
paste tmp.sample_name.txt PA-4_result.txt ${TMP_reference}/SGA_001/${sample_name}.stat.txt > ${sample_name}.SGA_001.result.txt
rm -f tmp.sample_name.txt

rm -f ${sample_name}.NanoFilt.BlastTab+.txt
