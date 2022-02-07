sample_name=`echo $1 | cut -d. -f1`
MODULE_HOME=$2
TMP_reference=${MODULE_HOME}/Reference/TMP/${sample_name}_analysis

mkdir ${sample_name}_SGA_003_analysis
cd ${sample_name}_SGA_003_analysis
mkdir -p ${TMP_reference}/SGA_003 #创建存放reference的临时文件

###################要分析点突变的基因##############################
cat ${MODULE_HOME}/Reference/SGA_003/* | dos2unix | seqkit seq -u > ${TMP_reference}/SGA_003/SGA_003.fasta #SGA_003下所有参考序列合并
java -jar /home/chensy/software/picard.jar CreateSequenceDictionary -REFERENCE ${TMP_reference}/SGA_003/SGA_003.fasta -OUTPUT ${TMP_reference}/SGA_003/SGA_003.dict #使用picard软件生成dict
lastdb ${TMP_reference}/SGA_003/SGA_003.db ${TMP_reference}/SGA_003/SGA_003.fasta #根据参考基因构建数据库
grep ">" ${TMP_reference}/SGA_003/SGA_003.fasta | sed "s/>//g" > ${TMP_reference}/SGA_003/reference_info_SGA_003.txt
##############################################################################

lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 ${TMP_reference}/SGA_003/SGA_003.db ../$1 | maf-convert -f ${TMP_reference}/SGA_003/SGA_003.dict sam -r 'ID:[id] PL:[nanopore] SM:[sample]' - > ${sample_name}.SGA_003.sam #last比对并且maf格式文件转换成sam格式文件
lastal -Q1 -P 4 -q 1 -b 1 -Q 0 -a 1 -e 45 -f BlastTab+ ${TMP_reference}/SGA_003/SGA_003.db ../$1 | grep -v ^# - > ${sample_name}.NanoFilt.BlastTab+ #比对生成blasttab+格式文件

grep "^@" ${sample_name}.SGA_003.sam > ${sample_name}.tmp1.sam #提取表头
grep -v "^@" ${sample_name}.SGA_003.sam > ${sample_name}.tmp2.sam #去除表头

if test -s ${sample_name}.tmp2.sam; then
Rscript ${MODULE_HOME}/sub/align_filter.R ${sample_name}.NanoFilt.BlastTab+ ${sample_name}.tmp2.sam ${sample_name}.NanoFilt.BlastTab+.txt ${sample_name}.tmp2.sam_filtered
else
touch ${sample_name}.NanoFilt.BlastTab+.txt
touch ${sample_name}.tmp2.sam_filtered
fi

cat ${sample_name}.tmp1.sam ${sample_name}.tmp2.sam_filtered > ${sample_name}.SGA_003.sam_filtered

rm -f ${sample_name}.tmp1.sam
rm -f ${sample_name}.tmp2.sam
rm -f ${sample_name}.tmp2.sam_filtered
rm -f ${sample_name}.SGA_003.sam
rm -f ${sample_name}.NanoFilt.BlastTab+
#rm -f ${sample_name}.NanoFilt.BlastTab+.txt

samtools view -bS ${sample_name}.SGA_003.sam_filtered | samtools sort - > ${sample_name}.SGA_003.sorted.bam #sam文件转化成bam文件并排序
samtools index ${sample_name}.SGA_003.sorted.bam #bam文件索引
samtools mpileup -a -s -f ${TMP_reference}/SGA_003/SGA_003.fasta ${sample_name}.SGA_003.sorted.bam > ${sample_name}.SGA_003.mpileup.txt

mkdir Fit_fasta
mkdir Mpileup_count

if test -s ${sample_name}.SGA_003.mpileup.txt; then
samtools mpileup -f ${TMP_reference}/SGA_003/SGA_003.fasta ${sample_name}.SGA_003.sorted.bam | java -jar ${MODULE_HOME}/sub/VarScan.v2.3.9.jar mpileup2indel --min-var-freq 0.5 –output-vcf 1 > ${sample_name}.SGA_003.InDel.vcf #Varscan生成插入缺失
python ${MODULE_HOME}/sub/SGA_003_mpileup.py -i ../$1 -d ${sample_name}.SGA_003.InDel.vcf
else
touch ${sample_name}.SGA_003.InDel.vcf
sed -i "1i Chrom\tPosition\tRef\tVar\tALT_frequency" ${sample_name}.SGA_003.InDel.vcf
fi

rm -f ${sample_name}.SGA_003.sam_filtered
rm -f ${sample_name}.SGA_003.sorted.bam
rm -f ${sample_name}.SGA_003.sorted.bam.bai
rm -f ${sample_name}.SGA_003.mpileup.txt


###统计测序深度

mkdir abstract_fastq
mkdir abstract_fasta

touch ${TMP_reference}/SGA_003/${sample_name}.stat.txt
echo "Depth" >> ${TMP_reference}/SGA_003/${sample_name}.stat.txt

echo "sample" > ${TMP_reference}/SGA_003/reads_stat.txt
cat ${TMP_reference}/SGA_003/reads_stat.txt ${TMP_reference}/SGA_003/reference_info_SGA_003.txt > ${TMP_reference}/SGA_003/reads_stat_info.txt

for i in $(cat ${TMP_reference}/SGA_003/reference_info_SGA_003.txt) #统计样本比对到的reads数
do
grep -w ${i} ${sample_name}.NanoFilt.BlastTab+.txt | awk '{print $1}' > ${sample_name}.${i}.abstract.txt #提取比对上的序列名

if test -s ${sample_name}.${i}.abstract.txt; then
seqkit grep -f ${sample_name}.${i}.abstract.txt ../$1 > ${sample_name}.${i}.abstract.fastq #提取比对上的序列
seqkit fq2fa ${sample_name}.${i}.abstract.fastq -o ${sample_name}.${i}.abstract.fasta #fastq转化成fasta
mv ${sample_name}.${i}.abstract.fastq abstract_fastq
mv ${sample_name}.${i}.abstract.fasta abstract_fasta

fi

let count=`cat ${sample_name}.${i}.abstract.txt | wc -l` #得到比对到的reads数
echo "$count" >> ${TMP_reference}/SGA_003/${sample_name}.stat.txt

rm -f ${sample_name}.${i}.abstract.txt
done

###和二级数据库比较
cat ${MODULE_HOME}/database/SGA_003/* | dos2unix | seqkit seq -u > ${TMP_reference}/SGA_003/SGA_003.database.fasta #SGA_003下所有参考序列合并
makeblastdb -in ${TMP_reference}/SGA_003/SGA_003.database.fasta -dbtype nucl -parse_seqids -hash_index -out ${TMP_reference}/SGA_003/SGA_003.database.db #建数据库


touch ${sample_name}.SGA_003_result.txt
echo "$sample_name" >> ${sample_name}.SGA_003_result.txt

for j in $(cat ${TMP_reference}/SGA_003/reference_info_SGA_003.txt) #统计样本比对到的reads数
do
if [ -f Fit_fasta/${sample_name}.${j}.fit.fasta ];then
blastn -query Fit_fasta/${sample_name}.${j}.fit.fasta -db ${TMP_reference}/SGA_003/SGA_003.database.db -evalue 1e-5 -outfmt 7 -num_alignments 1 -num_threads 2 | grep -v ^# - > ${j}.tmp.blast

if test -s ${j}.tmp.blast; then
awk '{print "best_matches("$2") identity("$3") mismatches("$5") gap("$6")"}' ${j}.tmp.blast >> ${sample_name}.SGA_003_result.txt
else
echo 'NM' >> ${sample_name}.SGA_003_result.txt
fi

else
echo '0' >> ${sample_name}.SGA_003_result.txt
fi
rm -f ${j}.tmp.blast
done

paste ${TMP_reference}/SGA_003/reads_stat_info.txt ${sample_name}.SGA_003_result.txt ${TMP_reference}/SGA_003/${sample_name}.stat.txt > ${sample_name}.SGA_003.result.txt
rm -f ${sample_name}.SGA_003_result.txt

rm -f ${sample_name}.NanoFilt.BlastTab+.txt
