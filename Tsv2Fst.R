#!/usr/bin/Rscript
library(fst)
library(data.table)
Tsv2Fst <- function(file){
  Output=sub(".tsv$",".fst",file)
  DT    <-fread(file,
                header = TRUE,
                sep='\t',
                col.names  = c("Reference","Position","Read"),
                colClasses = c("character","numeric","numeric"))
  write.fst(x=DT,path = Output,compress = 100)
}
INPUT='E:/Qsync/MY_script/GitHub/EJCAnalysis/Datasets/Arabidopsis/'
Cnt.file.cDNA = list.files(INPUT, pattern = "*.tsv", full.names = TRUE)
for (File in Cnt.file.cDNA){
  Tsv2Fst(File)
}
