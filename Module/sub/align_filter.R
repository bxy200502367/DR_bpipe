Args <- commandArgs()

# 加载包
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# 设置文件位置
blast_input_info <- Args[6]
sam_input_info <- Args[7]

#导入输入文件
blast_input <- fread(blast_input_info, header = FALSE, sep = "\t", stringsAsFactors = FALSE , check.names=FALSE)
sam_input <- fread(sam_input_info, header = FALSE, sep = "\t", stringsAsFactors = FALSE , check.names=FALSE)

#两者先排序
blast_input_sort <- blast_input[order(blast_input$V1,blast_input$V2),]
sam_input_sort <- sam_input[order(sam_input$V1,sam_input$V3),]

blast_input_sort$line_numbers <- c(1:nrow(blast_input_sort)) #增加一列，标记行号
blast_input_sort_filtered <- blast_input_sort %>% 
  dplyr::filter(V3>=85 & V4>=400) %>% 
  group_by(V1) %>% 
  dplyr::filter(V3 == max(V3))

sam_input_sort_filtered <- sam_input_sort[blast_input_sort_filtered$line_numbers,] #提取sam文件对应的行

write.table(blast_input_sort_filtered, Args[8], sep = "\t", col.names = FALSE ,row.names = FALSE, quote = FALSE)
write.table(sam_input_sort_filtered, Args[9], sep = "\t", col.names = FALSE ,row.names = FALSE, quote = FALSE)
