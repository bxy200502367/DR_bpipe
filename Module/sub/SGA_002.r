Args <- commandArgs()


# 加载包
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

blast_position_info <- Args[6]
blast_position_info_txt <- fread(blast_position_info, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_length_info <- Args[7]
gene_length_info_txt <- fread(gene_length_info, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

blast_length <- blast_position_info_txt$subject_end-blast_position_info_txt$subject_start+1
#table(blast_position_info_txt$subject_start)
df <- data.frame(Position = c(1:gene_length_info_txt$V2), 
                 Start_Frequency = NA,
				 End_Frequency = NA
				 )
df_tmp <- data.frame(Position = c(1:gene_length_info_txt$V2), 
                 Start_Frequency = NA,
				 End_Frequency = NA
				 )

df2 <- data.frame(Position = c(1:gene_length_info_txt$V2), 
                 length_Frequency = NA
				 )


for (i in 1:gene_length_info_txt$V2){
start1 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i-5),]
start2 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i-4),]
start3 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i-3),]
start4 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i-2),]
start5 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i-1),]
start6 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i),]
start7 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i+1),]
start8 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i+2),]
start9 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i+3),]
start10 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i+4),]
start11 <- blast_position_info_txt[which(blast_position_info_txt$subject_start == i+5),]
start_result <- rbind(start1,start2,start3,start4,start5,start6,start7,start8,start9,start10,start11)
start_result_length <- nrow(start_result)
start_position_frequency <- start_result_length/length(blast_length)
start_position_frequency_tmp <- round(start_position_frequency,4)*100
start_position_frequency <- paste0(round(start_position_frequency,4)*100,"%")

df[i,]$Start_Frequency <- start_position_frequency
df_tmp[i,]$Start_Frequency <- start_position_frequency_tmp

end1 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i-5),]
end2 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i-4),]
end3 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i-3),]
end4 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i-2),]
end5 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i-1),]
end6 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i),]
end7 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i+1),]
end8 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i+2),]
end9 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i+3),]
end10 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i+4),]
end11 <- blast_position_info_txt[which(blast_position_info_txt$subject_end == i+5),]
end_result <- rbind(end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,end11)
end_result_length <- nrow(end_result)
end_position_frequency <- end_result_length/length(blast_length)
end_position_frequency_tmp <- round(end_position_frequency,4)*100
end_position_frequency <- paste0(round(end_position_frequency,4)*100,"%")

df[i,]$End_Frequency <- end_position_frequency
df_tmp[i,]$End_Frequency <- end_position_frequency_tmp
}

write.table(df,Args[8], sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)

for (i in 1:gene_length_info_txt$V2){
length1 <- length(which(blast_length == i-5))
length2 <- length(which(blast_length == i-4))
length3 <- length(which(blast_length == i-3))
length4 <- length(which(blast_length == i-2))
length5 <- length(which(blast_length == i-1))
length6 <- length(which(blast_length == i))
length7 <- length(which(blast_length == i+1))
length8 <- length(which(blast_length == i+2))
length9 <- length(which(blast_length == i+3))
length10 <- length(which(blast_length == i+4))
length11 <- length(which(blast_length == i+5))
length_result <- length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11
length_frequency <- length_result/length(blast_length)
length_frequency <- paste0(round(length_frequency,4)*100,"%")
df2[i,]$length_Frequency <- length_frequency
}
write.table(df2, Args[9], sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)

start_judge <- mean(df_tmp$Start_Frequency[c(1:6)])
end_judge <- mean(df_tmp$End_Frequency[c((gene_length_info_txt$V2-5):gene_length_info_txt$V2)])
bread_point_position <- which.max(df_tmp$End_Frequency[c(7:(gene_length_info_txt$V2-6))])+6+4

if(start_judge+end_judge < 120){
shell_cmd <- paste0("echo 'M(",paste0(round(start_judge+end_judge,2),"%"),")' > tmp_wild_mutant.txt")
shell_cmd_2 <- paste0("echo 'pos",bread_point_position,"(",paste0(round(df_tmp$End_Frequency[bread_point_position],2),"%"),")'"," > tmp_break_point.txt")
system(shell_cmd,intern = TRUE)
system(shell_cmd_2,intern = TRUE)
}else{
shell_cmd <- paste0("echo 'W(",paste0(round(start_judge+end_judge,2),"%"),")' > tmp_wild_mutant.txt")
shell_cmd_2 <- paste0("echo '-' > tmp_break_point.txt")
system(shell_cmd,intern = TRUE)
system(shell_cmd_2,intern = TRUE)
}


