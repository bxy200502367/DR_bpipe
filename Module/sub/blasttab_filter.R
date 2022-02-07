Args <- commandArgs()

## 加载包
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

## 载入function
FilterBlasttab <- function(blasttabContent){
blasttabContent.filtered <- blasttabContent %>% 
  dplyr::filter(V3>=85 & V4>=400) %>% 
  group_by(V1) %>% 
  dplyr::filter(V3 == max(V3))
return(blasttabContent.filtered)
}

## 设置文件位置
blasttabInput <- Args[6]

##导入输入文件
blasttabContent <- fread(blasttabInput, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

blasttabContent.filtered <- FilterBlasttab(blasttabContent)

write.table(blasttabContent.filtered, Args[7], sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
