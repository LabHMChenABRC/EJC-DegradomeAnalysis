#!/usr/bin/Rscript
if(! require("data.table") | ! require("fst")){
  stop()
}
ARGS   = commandArgs(trailingOnly=TRUE)
if (length(ARGS)!=4) {
  stop("Four arguments must be supplied\nUseage:NMDTarget_EJC.R <Path_of_Dataset_folder> <GFF_file> <Representative_gene_model_file> <Output_of_Path>", call.=FALSE)
} else {
  INPUT   = ARGS[1]
  GFF_FILE= ARGS[2]
  REP_FILE= ARGS[3]
  OUTPUT  = ARGS[4]
  if (! file.exists(GFF_FILE)){
    stop("GFF file is required")
  }
  if (! file.exists(REP_FILE)){
    stop("Representative gene model file is required.")
  }
}
# test 
# INPUT   = "./Dataset/Arabidopsis/"
# GFF_FILE= "./Reference/Arabidopsis/TAIR10_GFF3_genes_transposons.gff"
# REP_FILE= "./Reference/Arabidopsis/Representative_gene_model_With_annotation.tsv"
# OUTPUT  = "./Output/Arabidopsis/"

# INPUT   = "./Dataset/Rice/"
# GFF_FILE= "./Reference/Rice/all.gff3"
# REP_FILE= "./Reference/Rice/Representative_gene_model_With_annotation.tsv"
# OUTPUT  = "./Output/Rice/"

FetchStackedCntFst <- function(files){
  CntTableList <- lapply(grep(".fst$",files,invert = FALSE,value = TRUE),
                         function(file){
                           SampleName=basename(sub("(.+)_.+.fst$","\\1",file))
                           read.fst(file,as.data.table=TRUE)[,"Sample":=SampleName
                                                             ][,.SD,.SDcols=c("Sample","Reference","Position","Read")]
                         })
  rbindlist(CntTableList)
}

Subset_MaxStackedcount<-function(DT){
  rbindlist(
    lapply(names(DT),function(x){
      DT[[x]][,.(Position=Position,Read=Read,Rank=rank(Read*-1, ties.method= "first")),by=c("Sample","Reference")][Rank==1][,`:=`(Type=x,Rank=NULL)][]
    })
  )
}

cDNA2Chr<-function(model,transcripts,positions){
  model_exon<-model[Transcript %in% transcripts][grep("exon",Type)][,exon_order:=-1][Type=="exon",exon_order:=as.numeric(1:.N),by=.(Transcript)]
  PosDT<-data.table(OriOrder=1:length(positions),"Transcript"=transcripts,"Position"=positions,"Position_dup"=positions)
  setkey(PosDT     ,Transcript,Position,Position_dup)
  setkey(model_exon,Transcript,Start_cDNA,End_cDNA)
  
  Hits<-foverlaps(PosDT,model_exon, type="within",nomatch=0L)[order(OriOrder)]
  Hits[Strand=="+",`:=`(Chr=Ref,Position_Chr=Position-Start_cDNA+Start)]
  Hits[Strand=="-",`:=`(Chr=Ref,Position_Chr=End_cDNA-Position+Start)]
  Hits[,.SD,.SDcols=c("Chr","Strand","Position_Chr","exon_order")]
  return(list("Chr"=Hits$Chr,"Strand"=Hits$Strand,"Position_Chr"=Hits$Position_Chr,"Gene"=Hits$Gene,"exon_order"=Hits$exon_order))
}

Chr2cDNA <- function(model,transcripts,positions){
  model_exon<-copy(model[Type=="exon"])
  PosDT<-data.table(OriOrder=1:length(positions),"Transcript"=transcripts,"Position"=positions,"Position_dup"=positions)
  setkey(PosDT     ,Transcript,Position,Position_dup)
  setkey(model_exon,Transcript,Start   ,End)
  # calculate model cDNA position
  Hits<-foverlaps(PosDT,model_exon, type="within",nomatch=0L)[order(OriOrder)]
  Hits[Strand=="+",Position_cDNA:=Position-Start+Start_cDNA]
  Hits[Strand=="-",Position_cDNA:=Start-Position+End_cDNA]
  return(Hits$Position_cDNA)
}

Add_UTR_Region <- function(model){
  Transcript_Len  <-model[Type == "exon"][,.(Len=sum(Len)),by=Transcript]
  setkey(Transcript_Len,Transcript)
  model_UTR5_Range<-model[Type == "CDS"
                          ][order(Transcript,Start_cDNA,End_cDNA),
                            ][,`:=`(last=.N,order=1:.N),by=.(Transcript)
                              ][order==1 & Start_cDNA > 1 
                                ][,`:=`(End_cDNA=Start_cDNA-1,
                                        Start_cDNA=1)
                                  ][,.SD,.SDcols=c("Transcript","Start_cDNA","End_cDNA")]
  
  model_UTR3_Range<-model[Type == "CDS"
                          ][order(Transcript,Start_cDNA,End_cDNA),
                            ][,`:=`(last=.N,order=1:.N),by=.(Transcript)
                              ][order==last
                                ][,`:=`(Start_cDNA=End_cDNA+1)
                                  ][,.SD,.SDcols=c("Transcript","Start_cDNA","End_cDNA")]
  model_UTR3_Range$End_cDNA<-Transcript_Len[.(model_UTR3_Range$Transcript),Len]
  model_UTR3_Range<-model_UTR3_Range[Start_cDNA<=End_cDNA]
  model_exon   <-model[Type == "exon"]
  
  setkey(model_exon,Transcript,Start_cDNA,End_cDNA)
  setkey(model_UTR5_Range,Transcript,Start_cDNA,End_cDNA)
  setkey(model_UTR3_Range,Transcript,Start_cDNA,End_cDNA)
  model_UTR5<-foverlaps(model_exon,model_UTR5_Range, type="any",nomatch=0L)
  model_UTR3<-foverlaps(model_exon,model_UTR3_Range, type="any",nomatch=0L)
  model_UTR<-rbindlist(
    list(
      model_UTR5[i.End_cDNA>End_cDNA,`:=`(i.End_cDNA=End_cDNA)
                 ][,`:=`(Start_cDNA=i.Start_cDNA,End_cDNA=i.End_cDNA,Type="UTR5")],
      
      model_UTR3[i.Start_cDNA<Start_cDNA,`:=`(i.Start_cDNA=Start_cDNA)
                 ][,`:=`(Start_cDNA=i.Start_cDNA,End_cDNA=i.End_cDNA,Type="UTR3")]
    )
  )[Strand=="+",`:=`(Start = cDNA2Chr(model,Transcript,Start_cDNA)$Position_Chr,
                     End   = cDNA2Chr(model,Transcript,End_cDNA)$Position_Chr)
    ][Strand=="-",`:=`(Start = cDNA2Chr(model,Transcript,End_cDNA)$Position_Chr,
                       End   = cDNA2Chr(model,Transcript,Start_cDNA)$Position_Chr)
      ][,Len := End_cDNA-Start_cDNA+1]
  model_merge<-rbindlist(
    list(model,
         model_UTR[,.SD,.SDcols=names(model)])
  )[,sort_model(.SD,.BY),by=.(Strand)
    ][order(Ref,Gene,Transcript)
      ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Len","Start_cDNA","End_cDNA")]
  return(model_merge)
}

Get_Attribute <- function(Type,Attributes,AttName){
  switch(
    Type,
    GTF = sub(sprintf('.*%s "(.*?)".*',AttName),'\\1',Attributes),
    GFF = sub(sprintf('.*%s=(.*?);.*',AttName),'\\1',Attributes)
  )
}

FiterTranscriptWithInCorrectCDS <- function(model){
  model_temp<-copy(model)
  setkey(model_temp,Transcript,Start,End)
  CDS              <- model_temp[Type == "CDS"]
  model_exon       <- model_temp[Type == "exon"]
  hit_row          <- foverlaps(CDS ,model_exon, type="within",nomatch=0L,which = TRUE)$xid
  CDSNotFullInExon <- CDS[!1:nrow(CDS) %in% hit_row,unique(Transcript)]
  return(model[! Transcript %in% CDSNotFullInExon])
}

Fetch_GFF3 <- function(model.file,Representative_Genes_models){
  GFF3   <-fread(model.file,header = FALSE,blank.lines.skip=TRUE,
                 colClasses = c("character","character","character","numeric","numeric","character","character","character","character"),
                 col.names = c("Ref","Source","Type","Start","End","Score","Strand","Phase","Attribute"))
  mRNA   <-GFF3[grep("mRNA|Transcript",Type,perl=TRUE)
                ][,{Temp<-unlist(strsplit(Attribute,split = ";"));
                .(Gene=sub("Parent=","",grep("Parent=",Temp,value=TRUE)),
                  Transcript=sub("ID=","",grep("ID=",Temp,value=TRUE)))}][]
  feature<-GFF3[grep("exon|CDS|five_prime_UTR|three_prime_UTR",Type,perl=TRUE,ignore.case = TRUE)
                ][grep("exon",Type,ignore.case = TRUE),`:=`(Type="exon")
                  ][grep("five_prime_UTR",Type,ignore.case = TRUE),Type:="UTR5"
                    ][grep("three_prime_UTR",Type,ignore.case = TRUE),Type:="UTR3"
                      ][,`:=`(Transcript=sub("Parent=","",grep("Parent=",unlist(strsplit(Attribute,split = ";")),value=TRUE)))]
  # Split feature when parent of feature contains multiple values
  if(feature[grep(",",Transcript),.N]){
    feature <- feature[, strsplit(Transcript, ",") ,by = names(feature)
                       ][,Transcript:=NULL
                         ][,setnames(.SD,"V1","Transcript")]  
  }
  
  GFF3   <-mRNA[feature,on="Transcript"
                ][,sort_model(.SD,.BY),by=.(Strand)
                  ][order(Ref,Gene,Transcript)
                    ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript")]
  # Protein-coding gene
  GFF3 <- GFF3[ Gene %in% GFF3[Type=="CDS",unique(Gene)]]
  # If Representative gene model contains one more intron.
  GFF3 <- GFF3[ Gene %in% GFF3[Transcript %in% Representative_Genes_models$Representative_gene_model
                            ][grep("exon",Type,ignore.case = TRUE)
                              ][,.(N=.N),by=.(Gene,Transcript)
                                ][N>1,unique(Gene)]]
  # Discard Gene contain CDS which is not within the exon region.(Some rice genes have mistakes)
  GFF3<-FiterTranscriptWithInCorrectCDS(GFF3)
  # Assign cDNA position
  GFF3.Exon<-GFF3[Type=="exon"
                ][,Len:=End-Start+1
                  ][,End_cDNA:=cumsum(Len),by=.(Gene,Transcript)
                    ][,Start_cDNA:=End_cDNA-Len+1
                      ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Len","Start_cDNA","End_cDNA")]
  GFF3.Features <-GFF3[Type!="exon"
                       ][,Len:=End-Start+1
                         ][Strand=="+",Start_cDNA:=Chr2cDNA(GFF3.Exon,Transcript,Start)
                           ][Strand=="-",Start_cDNA:=Chr2cDNA(GFF3.Exon,Transcript,End)
                             ][,End_cDNA:=Start_cDNA+Len-1]
  model <- rbindlist(list(GFF3.Exon,GFF3.Features))
  
  if(model[Type %in% c("UTR5","UTR3"),.N]==0){
    model <- Add_UTR_Region(model)
  }
  model<-model[order(Ref,Gene,Transcript,Start_cDNA,End_cDNA),
        ][,isoform_num:=length(unique(Transcript)),by=.(Gene)
          ][,`:=`(num=.N,order=1:.N),by=.(Transcript,Type)
            ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Start_cDNA","End_cDNA","Len","isoform_num","num","order")]
}

EJC_Region <- function(model,from=30,to=26){
  if ( from < to ){
    Stop("from should >= to")
  }else{
    copy(model)[Type=="exon" & order<num & Len>to
                ][,`:=`(Start_cDNA = End_cDNA - from,
                        End_cDNA   = End_cDNA - to,
                        Type       ="EJC")
                  ][Start_cDNA<=0,Start_cDNA:=1
                    ][Strand=="+",`:=`(Start=cDNA2Chr(model,Transcript,Start_cDNA)$Position_Chr,
                                       End  =cDNA2Chr(model,Transcript,End_cDNA)$Position_Chr)
                      ][Strand=="-",`:=`(Start=cDNA2Chr(model,Transcript,End_cDNA)$Position_Chr,
                                         End  =cDNA2Chr(model,Transcript,Start_cDNA)$Position_Chr)
                        ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Start_cDNA","End_cDNA")]
  }
}

PeakinEJC<-function(DT,Feature_model,CDS,UTR){
  DTname=c(names(DT))
  DT<-copy(DT)[,Position_dup:=Position]
  setkey(DT,Reference,Position,Position_dup)
  setkey(Feature_model,Transcript,Start_cDNA,End_cDNA)
  hit_index<-foverlaps(DT,Feature_model,type = "within",nomatch=0L)[,.SD,.SDcols=c(DTname,"Start","End","Position_dup")
                                                                    ][,setnames(.SD,c("Start","End"),c("Start_EJC","End_EJC"))
                                                                      ][,c("Chr","Strand","Position_Chr","Gene","order"):=cDNA2Chr(Model,Reference,Position)
                                                                        ][order(Reference,Position),.SD,by=.(Sample,Gene)]
  Model_cols=c("Transcript","Type","Start_cDNA","End_cDNA")
  FullModel<-rbindlist(list(CDS[,.SD,.SDcols=Model_cols],
                            UTR[,.SD,.SDcols=Model_cols]))
  setkey(FullModel,Transcript,Start_cDNA,End_cDNA)
  setkey(hit_index,Reference ,Position,Position_dup)
  hit_index_ann<-foverlaps(hit_index,FullModel,type = "within",nomatch=0L)[,.SD,.SDcols=c(DTname[1],"Gene",DTname[2:length(DTname)],"Position_Chr")]
  return(hit_index_ann)
}

SubsetByFeature<-function(model,DTmapping){
  DT<-DTmapping[Reference %in% unique(model$Transcript)]
  DT[,Position_dup:=Position]
  setkey(DT,Reference,Position,Position_dup)
  setkey(model,Transcript,Start_cDNA,End_cDNA)
  Hits<-foverlaps(DT,model, type="within",nomatch=0L)
  return(Hits[,.SD,.SDcols=names(DTmapping)])
}

Build_AS_Model<-function(CntDT,model,CDS){
  Gene_exp_DT<- CntDT[,unique(.SD),.SDcols=c("Gene","Reference")][,N:=1:.N][]
  Model_exp  <- copy(model)[Gene %in% Gene_exp_DT$Gene]
  CDS_exp    <- copy(CDS)[Gene %in% Gene_exp_DT$Gene]
  setkey(Model_exp,Gene)
  setkey(Gene_exp_DT,Gene)
  setkey(CDS_exp,Gene)
  AS.DT<-rbindlist(lapply(Gene_exp_DT$Gene,
                          function(x) AnnotatedAS(Model_exp[.(x)],CDS_exp[.(x)],Gene_exp_DT[.(x),Reference]) ))
  Output_colnames=names(AS.DT)
  rbindlist(list(
    AS.DT,
    Model_exp[Transcript %in% Gene_exp_DT$Reference
              ][,AS:=0][,.SD,.SDcols=Output_colnames]  
  ))[,sort_model(.SD,.BY),by=.(Strand)
     ][order(Ref,Gene,Transcript)
       ][,.SD,.SDcols=Output_colnames]
  
}

AnnotatedAS<-function(model,CDS,Major){
  Model_Majorform <- copy(model[Transcript == Major,])
  Output_colnames  = c("Ref","Type","Start","End","Strand","Gene","Transcript","AS")
  # Major Gene form should contain intron and total isoforms of Gene > 1
  if(Model_Majorform$isoform_num==1 || Model_Majorform$num==1){
    return(
      setNames(data.table(matrix(nrow=0,ncol=8)),Output_colnames)
    )
  }
  Model_Isoforms      <- copy(model[!Transcript %in% Major,])
  CDS_Majorform_frist <- CDS[Transcript == Major & order==1,]
  CDS_Isoforms_frist  <- CDS[!Transcript %in% Major & order==1 ,]
  CDS_Majorform_last  <- CDS[Transcript == Major & order==num,]
  CDS_Isoforms_last   <- CDS[!Transcript %in% Major & order==num ,]
  if (Model_Majorform[1,Strand]=="+") {
    Major_start = CDS_Majorform_frist[,Start]
    Major_end   = CDS_Majorform_last[,End]
    CDS_dif     = unique(c(CDS_Isoforms_frist[! Start %in% Major_start,Transcript],
                           CDS_Isoforms_last[! End   %in% Major_end  ,Transcript]))
  } else {
    Major_start = CDS_Majorform_frist[,End]
    Major_end   = CDS_Majorform_last[,Start]
    CDS_dif     = unique(c(CDS_Isoforms_frist[! End  %in% Major_start,Transcript],
                           CDS_Isoforms_last[! Start %in% Major_end  ,Transcript]))
  }
  Model_Isoforms<-Model_Isoforms[ Transcript %in% CDS_dif]
  
  AS_ann<-lapply(split(Model_Isoforms,by="Transcript"),function(x){
    AS<-rbindlist(list(Model_Majorform,x))[ End>=max(x[,min(Start)],Model_Majorform[,min(Start)]) & 
                                              Start<=min(x[,max(End)]  ,Model_Majorform[,max(End)]) 
                                            ][,`:=`(num=.N,order=1:.N),by=.(Transcript)][]
    if (Model_Majorform[1,Strand]=="+") {
      AS<-AS[ ,`:=`(intron_3ss=.N),by=.(Gene,Start)
              ][ ,`:=`(intron_5ss=.N),by=.(Gene,End) 
                 ][ order==1,intron_3ss:=2
                    ][ order==num,intron_5ss:=2][]
      AS_ann<-rbindlist(list(
        AS[ intron_3ss==1 & intron_5ss==2 ][,`:=`(End  =Start,Type="3ss")],
        AS[ intron_3ss==2 & intron_5ss==1 ][,`:=`(Start=End  ,Type="5ss")],
        AS[ intron_3ss + intron_5ss==2 ][,`:=`(Type=ifelse(Transcript==Major,"Exon_Gain","Exon_Loss"))]
      ))[,`:=`(Transcript=Major,AS=1)][,.SD,.SDcols=Output_colnames]
    } else {
      AS<-AS[ ,`:=`(intron_3ss=.N),by=.(Gene,End)
              ][ , `:=`(intron_5ss=.N),by=.(Gene,Start)
                 ][ order==1,intron_3ss:=2
                    ][ order==num,intron_5ss:=2][]
      AS_ann<-rbindlist(list(
        AS[ intron_3ss==1 & intron_5ss==2 ][,`:=`(Start=End   ,Type="3ss")],
        AS[ intron_3ss==2 & intron_5ss==1 ][,`:=`(End  =Start ,Type="5ss")],
        AS[ intron_3ss + intron_5ss==2 ][,`:=`(Type=ifelse(Transcript==Major,"Exon_Gain","Exon_Loss"))]
      ))[,`:=`(Transcript=Major,AS=1)][,.SD,.SDcols=Output_colnames]
    }
  })
  return(rbindlist(AS_ann))
}

Build_STOP_UTR_intron_Model<-function(CntDT,model){
  Output_colnames  <- c("Ref","Type","Start","End","Strand","Gene","Transcript","STOP")
  Gene_exp_DT      <- CntDT[Type=="UTR3"][,unique(.SD),.SDcols=c("Reference")]
  model            <- copy(model[ Transcript %in% Gene_exp_DT$Reference ])
  STOP             <- rbindlist(
                       list(model[Strand=="+" & Type=="CDS" & num==order][,`:=`(Type="STOP",Start=End)],
                            model[Strand=="-" & Type=="CDS" & num==order][,`:=`(Type="STOP",End=Start)]))
  model_STOP  <- rbindlist(
                  list(model[Type=="exon"][,STOP:=0],
                       STOP[,STOP:=1])
                  )[,sort_model(.SD,.BY),by=.(Strand)
                    ][order(Ref,Gene,Transcript)
                      ][,.SD,.SDcols=Output_colnames]
  
}

sort_model<-function(x,strand){
  if(strand=="+"){
    return(x[order(End,Start,decreasing=FALSE),])
  }else{
    return(x[order(Start,End,decreasing=TRUE),] )
  }
}

NMD_check<-function(NMD,Pos,DistanceAllow=2){
  if(NMD[,max(AS)]==0){
    return(FALSE)
  }
  
  if(NMD[1,Strand=="+"]){
    temp <- NMD[Start<Pos] 
  }else{
    temp <- NMD[End>Pos]
  }
  
  N=0
  if ( nrow(temp)>0){
    if(max(temp$AS)>0){
      N=temp[,AS_cumsum:=rev(cumsum(rev(AS)))][AS_cumsum==0,.N]
    }
  }
  ifelse( N >= 1 & N <= DistanceAllow,
          TRUE,
          FALSE)
}

NMD_STOP_check<-function(Model.STOP,Pos,DistanceAllow=2,DistanceSTOP=50){
  if(Model.STOP[,max(STOP)]==0){
    return(FALSE)
  }
  if(Model.STOP[1,Strand=="+"]){
    temp <- Model.STOP[Start<Pos] 
  }else{
    temp <- Model.STOP[End>Pos]
  }
  Downstream_N=0
  if ( nrow(temp)>0){
    if(max(temp$STOP)>0){
      Downstream=temp[,STOP_cumsum:=rev(cumsum(rev(STOP)))][STOP_cumsum==0,]
      Downstream_N=nrow(Downstream)
    }
  }
  if( Downstream_N >= 1 & Downstream_N <= DistanceAllow){
    STOP_Pos=Model.STOP[Type=="STOP",Start]
    Exon_End_Pos=ifelse( Downstream[Downstream_N,Strand=="+"],
                         Downstream[Downstream_N,End],
                         Downstream[Downstream_N,Start])
    ifelse( abs(Exon_End_Pos-STOP_Pos) >= DistanceSTOP,
            TRUE,
            FALSE )
  }else{
    FALSE
  }
}

Summary_NMD_result<-function(Count_AS_DT){
  DT<-Count_AS_DT[,.SD,.SDcols=c("Gene","Reference","AS","UTR3","Read")
                  ][,Read:=max(Read),by=.(Gene)
                    ][,unique(.SD)
                      ][,NMD:=ifelse(AS==TRUE & UTR3==TRUE,"AS/Intron+ 3′UTR",
                                     ifelse(AS==TRUE,"AS","Intron+ 3′UTR"))
                        ][,.SD,.SDcols=c("Gene","Reference","NMD","Read")
                          ][,setnames(.SD,c("Locus","Representative gene model","Putative NMD trigger","Abundance of Max. 5′P peak (TP40M)"))]
}
######## Main ########
Max_Read_Cutoff = 10 # Max peak within transcript should >=10 TP40M
if(!dir.exists(OUTPUT)){
  dir.create(OUTPUT,recursive = TRUE)
}

# Annotation file
Gene_annotation            <-fread(REP_FILE,header=TRUE,sep="\t",drop = c("Representative_gene_model"))
Representative_Genes_models<-fread(REP_FILE,header=TRUE,sep="\t",select = c("Locus","Representative_gene_model"))

# Gene model
# Representative Genes models is intron-contain and protein-coding 
Model      <- Fetch_GFF3(GFF_FILE,Representative_Genes_models)
Model_EJC  <- EJC_Region(Model,from=30,to=26)
Model_exon <- Model[Type == "exon"]
Model_CDS  <- Model[Type == "CDS"]
Model_UTR  <- Model[Type %in% c("UTR5","UTR3")]

# Count table
Cnt.file.cDNA = list.files(INPUT, pattern = "*.fst",full.names = TRUE)
Stacked_count = list()
# Peaks from exon of representative genes models only
Stacked_count$exon    <-FetchStackedCntFst(Cnt.file.cDNA)[,Position:=Position-100 # move frist position to TSS
                                                          ][Reference %in% Representative_Genes_models$Representative_gene_model
                                                            ][,SubsetByFeature(Model_exon,.SD)]

# Max peak within cDNA
Stacked_count_max     <- Subset_MaxStackedcount(Stacked_count)
# Max peak within cDNA is located EJC region
Stacked_count_max_EJC <- PeakinEJC(Stacked_count_max,Model_EJC,Model_CDS,Model_UTR)

##### frist case(AS): representative genes model has AS events upstream of max peak(>=10 TP40M). ######
# Create Representative_Genes_models with/without AS annotation
Model_AS     <- Build_AS_Model(Stacked_count_max_EJC[Read>=eval(Max_Read_Cutoff)],Model_exon,Model_CDS)
setkey(Model_AS,Transcript)

##### second case(Intron+ 3′UTR): representative genes model has normal STOP >= 50-nt upstream of exon which has max peak(>=10 TP40M). ######
# Create Representative_Genes_models with STOP annotation
Model_STOP   <- Build_STOP_UTR_intron_Model(Stacked_count_max_EJC[Read>=eval(Max_Read_Cutoff)],Model)
setkey(Model_STOP,Transcript)

# Check those gene has AS events/STOP codon upstream of EJC peak, and there is <= one intron between EJC peak and AS events/STOP codon.
Stacked_count_max_EJC_AS<-
  copy(Stacked_count_max_EJC)[,`:=`(AS_event=ifelse(Read>=eval(Max_Read_Cutoff), # max peak >= 10 TP40M
                                                      NMD_check(Model_AS[.(Reference)],Position_Chr,2), # check AS
                                                      FALSE),
                                      UTR3_event=FALSE),
                                by=.(Sample,Reference,Position_Chr)
                                ][Type=="UTR3",
                                  UTR3_event:=ifelse(Read>=eval(Max_Read_Cutoff), # max peak >= 10 TP40M
                                                     NMD_STOP_check(Model_STOP[.(Reference)],Position_Chr,2,50), # check STOP
                                                     FALSE),
                                  by=.(Sample,Reference,Position_Chr)
                                  ][,`:=`(AS  =any(AS_event),
                                          UTR3=any(UTR3_event),
                                          NMD =any(AS_event,UTR3_event) ),
                                      by=.(Reference)
                                      ][NMD==TRUE,]

####  Output NMD Genes list
Fin          <- Summary_NMD_result(Stacked_count_max_EJC_AS)
# Merge Annonation
if( any(names(Gene_annotation) %in% "Q_val_up") ){
  Gene_annotation<-Gene_annotation[,Q_val_up:=as.character(Q_val_up)][is.na(Q_val_up),Q_val_up:="NA"][,setnames(.SD,"Q_val_up","Q val up(lba1 upf3-1/WT)")]
}
FinAnn       <-Gene_annotation[Fin,on="Locus"][,.SD,by=c("Locus","Representative gene model")]
# Ouput
fwrite(file = paste0(c(OUTPUT,"NMD_targets.tsv"),collapse = "/"),FinAnn,sep="\t")
####  Summary of NMD Genes based on AS/Intron+ 3′UTR
# Summay_NMD   <- Stacked_count_max_EJC_AS[,.(Gene_number=length(unique(Gene))),by=.(AS,UTR3)]
# fwrite(file = paste0(c(OUTPUT,"NMD_targets_Summay.tsv"),collapse = "/"),Summay_NMD,sep="\t")
