sample_name=`echo $1 | cut -d. -f1`
MODULE_HOME=$2
TMP_reference=${MODULE_HOME}/Reference/TMP/${sample_name}_analysis

mkdir ${sample_name}_SGA_002_analysis
cd ${sample_name}_SGA_002_analysis
mkdir -p ${TMP_reference}/SGA_002 #创建存放reference的临时文件

seqkit fq2fa ../$1 -o ${sample_name}.fasta #使用seqkit软件把fastq格式转换成fasta格式
cat ${MODULE_HOME}/Reference/SGA_002/* | dos2unix | seqkit seq -u > ${TMP_reference}/SGA_002/SGA_002.fasta #reference文件夹下的所有参考序列合并
makeblastdb -in ${TMP_reference}/SGA_002/SGA_002.fasta -dbtype nucl -parse_seqids -out ${TMP_reference}/SGA_002/SGA_002 #对adeNf参考基因建库
blastn -num_threads 12 -query ${sample_name}.fasta -db ${TMP_reference}/SGA_002/SGA_002 -out SGA_002_${sample_name}.out -outfmt '7 qacc sacc qstart qend sstart send' -evalue 1e-5
grep -v "#" SGA_002_${sample_name}.out | awk '{if ($5 > $6) {$7=$6;$8=$5} else {$7=$5;$8=$6};printf "%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$7,$8}' | sed '1i query_id\tsubject_id\tquery_start\tquery_end\tsubject_start\tsubject_end' > ${sample_name}.SGA_002.info #排序

rm -f SGA_002_${sample_name}.out

grep ">" ${TMP_reference}/SGA_002/SGA_002.fasta | sed "s/>//g" > ${TMP_reference}/SGA_002/reference_info_SGA_002.txt
seqkit fx2tab -l -g -n -i -H ${TMP_reference}/SGA_002/SGA_002.fasta > SGA_002_gene_length.txt ###统计reference每个基因长度

for i in $(cat ${TMP_reference}/SGA_002/reference_info_SGA_002.txt) #统计样本比对到的reads数
do
grep -w ${i} ${sample_name}.SGA_002.info | sed '1i query_id\tsubject_id\tquery_start\tquery_end\tsubject_start\tsubject_end' > ${sample_name}.${i}.SGA_002.info #提取比对上的序列
grep -w ${i} SGA_002_gene_length.txt > SGA_002.${i}.length_info

##生成各项指标文件
touch ${sample_name}.${i}.SGA_002_sample.txt
touch ${sample_name}.${i}.SGA_002_depth.txt
touch ${sample_name}.${i}.SGA_002_wild_mutant.txt
touch ${sample_name}.${i}.SGA_002_break_point.txt
echo "$sample_name" >> ${sample_name}.${i}.SGA_002_sample.txt ###样本名
let count=`cat ${sample_name}.${i}.SGA_002.info | wc -l`
echo "$count" >> ${sample_name}.${i}.SGA_002_depth.txt ###深度
if [ $count -ge 10 ] ; then
Rscript ${MODULE_HOME}/sub/SGA_002.r ${sample_name}.${i}.SGA_002.info SGA_002.${i}.length_info ${sample_name}.${i}.position_info ${sample_name}.${i}.length_info
cat tmp_wild_mutant.txt >> ${sample_name}.${i}.SGA_002_wild_mutant.txt
cat tmp_break_point.txt >> ${sample_name}.${i}.SGA_002_break_point.txt
rm -f tmp_wild_mutant.txt
rm -f tmp_break_point.txt
else
echo "LD" >> ${sample_name}.${i}.SGA_002_wild_mutant.txt
echo "-" >> ${sample_name}.${i}.SGA_002_break_point.txt
fi

##合并结果
paste ${sample_name}.${i}.SGA_002_sample.txt ${sample_name}.${i}.SGA_002_depth.txt ${sample_name}.${i}.SGA_002_wild_mutant.txt ${sample_name}.${i}.SGA_002_break_point.txt > ${sample_name}.${i}_result.txt

sed -i '1i '$i'\t'$i'_depth\t'$i'_result\tbreak_point' ${sample_name}.${i}_result.txt

rm -f ${sample_name}.${i}.SGA_002_sample.txt
rm -f ${sample_name}.${i}.SGA_002_depth.txt
rm -f ${sample_name}.${i}.SGA_002_wild_mutant.txt
rm -f ${sample_name}.${i}.SGA_002_break_point.txt
rm -f SGA_002.${i}.length_info
rm -f ${sample_name}.${i}.SGA_002.info
done

rm -f ${sample_name}.fasta
rm -f ${sample_name}.SGA_002.info
rm -f SGA_002_gene_length.txt

paste *_result.txt > ${sample_name}.SGA_002_result.txt
