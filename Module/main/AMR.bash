sample_name=`echo $1 | cut -d. -f1` #获取样本名
MODULE_HOME=$2
TMP_reference=${MODULE_HOME}/Reference/TMP/${sample_name}_analysis

mkdir ${sample_name}_AMR_analysis
cd ${sample_name}_AMR_analysis
mkdir -p ${TMP_reference}/AMR #创建存放reference的临时文件

##处理和整理reference
cat ${MODULE_HOME}/Reference/AMR/* | dos2unix | seqkit seq -u > ${TMP_reference}/AMR/AMR.fasta #reference文件夹下的所有参考序列合并
lastdb ${TMP_reference}/AMR/AMR.db ${TMP_reference}/AMR/AMR.fasta #根据参考基因构建数据库

lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 -f BlastTab+ ${TMP_reference}/AMR/AMR.db ../$1 | grep -v ^# - > ${sample_name}.NanoFilt.BlastTab+ #比对生成blasttab+格式文件
Rscript ${MODULE_HOME}/sub/blasttab_filter.R ${sample_name}.NanoFilt.BlastTab+ ${sample_name}.NanoFilt.BlastTab+.txt #过滤脚本过滤比对到的数据
rm -f ${sample_name}.NanoFilt.BlastTab+


##生成各项指标文件
touch ${TMP_reference}/AMR/${sample_name}.stat.txt
touch ${TMP_reference}/AMR/${sample_name}.stat_length.txt
touch ${TMP_reference}/AMR/${sample_name}.stat_score.txt
touch ${TMP_reference}/AMR/${sample_name}.stat_identity.txt
echo "$sample_name" >> ${TMP_reference}/AMR/$sample_name.stat.txt
echo "Ave_Length" >> ${TMP_reference}/AMR/$sample_name.stat_length.txt
echo "Ave_Score" >> ${TMP_reference}/AMR/$sample_name.stat_score.txt
echo "Ave_Identity" >> ${TMP_reference}/AMR/$sample_name.stat_identity.txt

echo "sample" > ${TMP_reference}/AMR/reads_stat.txt
grep ">" ${TMP_reference}/AMR/AMR.fasta | sed "s/>//g" > ${TMP_reference}/AMR/reference_info_AMR.txt
cat ${TMP_reference}/AMR/reads_stat.txt ${TMP_reference}/AMR/reference_info_AMR.txt > ${TMP_reference}/AMR/reads_stat_info.txt

for i in $(cat ${TMP_reference}/AMR/reference_info_AMR.txt) #统计样本比对到的reads数
do
grep -w ${i} ${sample_name}.NanoFilt.BlastTab+.txt | awk '{print $1}' > ${sample_name}.${i}.abstract.txt #提取比对上的序列名

if test -s ${sample_name}.${i}.abstract.txt; then
seqkit grep -f ${sample_name}.${i}.abstract.txt ../$1 > ${sample_name}.${i}.abstract.fastq #提取比对上的序列
seqkit fq2fa ${sample_name}.${i}.abstract.fastq -o ${sample_name}.${i}.abstract.fasta #fastq转化成fasta
else
touch ${sample_name}.${i}.abstract.fastq
touch ${sample_name}.${i}.abstract.fasta
fi
seqkit fx2tab -H -n -i -l -q ${sample_name}.${i}.abstract.fasta > ${sample_name}.${i}_length.xls #得到每一条测序read的长度信息
echo `awk '{print $2}' ${sample_name}.${i}_length.xls | awk '{sum+=$1} END {if(sum != 0){print sum/(NR-1)}}'` >> ${TMP_reference}/AMR/${sample_name}.stat_length.txt #得到平均测序长度
grep -w ${i} ${sample_name}.NanoFilt.BlastTab+.txt | awk '{print $1,$3,$15,$13}' > ${i}.${sample_name}.tmp.txt
echo `awk '{print $3}' ${i}.${sample_name}.tmp.txt | awk '{sum+=$1} END {if(sum != 0){print sum/NR}}'` >> ${TMP_reference}/AMR/${sample_name}.stat_score.txt #得到平均比对分值
echo `awk '{print $2}' ${i}.${sample_name}.tmp.txt | awk '{sum+=$1} END {if(sum != 0){print sum/NR}}'` >> ${TMP_reference}/AMR/${sample_name}.stat_identity.txt #得到平均比对identity值
let count=`cat ${i}.${sample_name}.tmp.txt | wc -l` #得到比对到的reads数
echo "$count" >> ${TMP_reference}/AMR/${sample_name}.stat.txt

rm -f ${sample_name}.${i}.abstract.txt
rm -f ${sample_name}.${i}_length.xls
rm -f ${i}.${sample_name}.tmp.txt

done


mkdir abstract_fastq
mkdir abstract_fasta
mv *abstract.fastq abstract_fastq
mv *abstract.fasta abstract_fasta

##合并结果
paste ${TMP_reference}/AMR/reads_stat_info.txt ${TMP_reference}/AMR/${sample_name}.stat.txt ${TMP_reference}/AMR/${sample_name}.stat_length.txt ${TMP_reference}/AMR/${sample_name}.stat_score.txt ${TMP_reference}/AMR/${sample_name}.stat_identity.txt > ${sample_name}.AMR.result_detail.txt
paste ${TMP_reference}/AMR/reads_stat_info.txt ${TMP_reference}/AMR/${sample_name}.stat.txt > ${sample_name}.AMR.result.txt

rm -f ${sample_name}.NanoFilt.BlastTab+.txt
