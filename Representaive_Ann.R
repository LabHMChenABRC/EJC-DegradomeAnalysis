library(data.table)

INPUT       = "."
OUTPUT_Rice = paste0(INPUT,"/Reference/Rice/Representative_gene_model_With_annotation.tsv")
OUTPUT_Ath  = paste0(INPUT,"/Reference/Arabidopsis/Representative_gene_model_With_annotation.tsv")

setwd(INPUT)
#### Rice Representative_gene_model_With_annotation ####
# Reference annotation is from MSU Release 7 database
Representative_gene_file_Rice   <-paste0(INPUT,"/Reference/Rice/all.locus_brief_info.7.0")
Gene_annotation_Rice            <-fread(Representative_gene_file_Rice,header=TRUE,sep="\t")[is_representative=="Y"
                                                                                            ][,setnames(.SD,
                                                                                                        c("locus","model","annotation"),
                                                                                                        c("Locus","Representative_gene_model","Description"))
                                                                                              ][,.SD,.SDcols=c("Locus","Representative_gene_model","Description")]
fwrite(x = Gene_annotation_Rice,file = OUTPUT_Rice,sep = "\t")

#### Arabidopsis Representative_gene_model_With_annotation ####

# Reference annotation from TAIR database
Gene_symbol_full_name_file <-paste0(INPUT,"/Reference/Arabidopsis/gene_aliases_20180930.txt")
Gene_function_file         <-paste0(INPUT,"/Reference/Arabidopsis/Araport11_functional_descriptions_20180330.txt")
Representative_gene_file   <-paste0(INPUT,"/Reference/Arabidopsis/TAIR10_representative_gene_models.txt")

# AS Dateset is from "(B)WT_vs_lba1upf3-1" sheet of Supplemental Data Set 1 in Drechsel et al. (2013). Plant Cell 10.1105/tpc.113.115485.
AS_UP_file                 <-paste0(INPUT,"/Reference/Arabidopsis/tpc115485Supp-WT_vs_lba1upf3-1.tsv")

Gene_function              <-fread(Gene_function_file,header=TRUE,sep="\t")[gene_model_type=="protein_coding",
                                                                            ][,setnames(.SD,c("name","Computational_description"),c("Locus","Description"))
                                                                              ][grep("^AT.*G.*",Locus)
                                                                                ][,Locus:=sub("(^AT.*G.*).\\d+","\\1",Locus)
                                                                                  ][,.SD[1],by=.(Locus)
                                                                                    ][,Description:=gsub(";*\\(source:.+\\)|;*\\(Source:.+\\)| *$","",Description,perl = TRUE)]

Representative_Genes_models<-fread(Representative_gene_file,header=FALSE,sep="\t",skip = 1)[,setnames(.SD,"Representative_gene_model")
                                                                                            ][,Locus:=sub("(^AT.*G.*).\\d+","\\1",Representative_gene_model)
                                                                                              ][Locus %in% Gene_function$Locus,
                                                                                                ][,.SD,.SDcols=c("Locus","Representative_gene_model")]
Gene_symbol                <-fread(Gene_symbol_full_name_file,header=TRUE,sep="\t")[name %in% Gene_function$Locus,
                                                                                    ][,.SD[1],by=.(name)
                                                                                      ][,setnames(.SD,c("name","symbol"),c("Locus","Symbol"))
                                                                                        ][grep("^AT.*G.*",Locus)]
AltSig                     <-fread(AS_UP_file,sep="\t",header=TRUE)[,setnames(.SD,names(.SD),gsub(" ","_",names(.SD)))
                                                                    ][,.(Q_val_up=min(Q_val_up)),by=.(gene_name)
                                                                      ][,setnames(.SD,"gene_name","Locus")]
Gene_annotation            <-Gene_symbol[Representative_Genes_models,on="Locus"]
Gene_annotation            <-Gene_function[Gene_annotation,on="Locus"]
Gene_annotation_Ath        <-AltSig[Gene_annotation,on="Locus"
                                    ][,Q_val_up:=as.character(Q_val_up)
                                      ][is.na(Q_val_up),Q_val_up:="NA"
                                        ][is.na(Symbol) | is.null(Symbol) ,Symbol:=""
                                          ][is.na(Description)|Description=="None",Description:=gene_model_type
                                            ][,.SD,.SDcols=c("Locus","Representative_gene_model","Symbol","Q_val_up","Description")]

fwrite(x = Gene_annotation_Ath,file = OUTPUT_Ath,sep = "\t")

