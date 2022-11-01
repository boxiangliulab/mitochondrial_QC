
setwd("D:/Cell Paper/Original image/Figure 1/Figure 1 I/two portion Z test/two portion Z test")
rm(list=ls())

library(readxl)
library(dplyr)
library(tidyr)

filename <- "PIVI-mitochondrial Area sizes 2.xlsx"
sheet_names = excel_sheets(path = filename)

sheet_name=sheet_names[2]
my_data <- read_excel(filename,sheet = sheet_name)

my_list=list()
for(i in 1:length(colnames(my_data))){
  nam <- paste("Group_", LETTERS[i], sep = "")
  assign(nam, as.vector(pull(drop_na(my_data[i]))))
  my_list[[i]] = get(nam)
}

Z_test=c()
control =  my_list[[1]]

for (i in 2:length(names(my_data))){
  group = my_list[[i]]
  for (threshold in 2:10){
    p1 = mean(control > threshold)
    p2 = mean(group  > threshold)
    all_data = c(control,group)
    p0 =mean(all_data  > threshold)
    Z_stat = (p1-p2)/sqrt(p0*(1-p0)*(1/length(control)+1/length(group)))
    temp = pnorm(abs(Z_stat),lower.tail = FALSE)*2
    Z_test =c(Z_test,(temp))}
}


PVALUESZ= matrix(Z_test,ncol = length(names(my_data))-1 ,byrow = FALSE)
writeClipboard(as.character((PVALUESZ)),format = 1)
write.table(PVALUESZ,col.names = names(my_data)[seq(2,length(names(my_data)))]
,"results.csv",sep=",", row.names = F)

