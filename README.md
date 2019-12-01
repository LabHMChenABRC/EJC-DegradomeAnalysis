# EJC-DegradomeAnalysis
 
**The putaive NMD target prediction :**

Execute R script "NMDTarget_EJC.R" will output putative NMD targets list, and it requires R package data.table and fst.
```
Rscript NMDTarget_EJC.R <Path_of_Dataset_folder> <GFF_file> <Representative_gene_model_file> <Output_of_Path>
```
```
#Supplemental Dataset 1 :
Rscript NMDTarget_EJC.R ./Dataset/Arabidopsis/ ./Reference/Arabidopsis/TAIR10_GFF3_genes_transposons.gff ./Reference/Arabidopsis/Representative_gene_model_With_annotation.tsv ./Ouput/Arabidopsis/


#Supplemental Dataset 2 :
Rscript NMDTarget_EJC.R ./Dataset/Rice/ ./Reference/Rice/all.gff3 ./Reference/Rice/Representative_gene_model_With_annotation.tsv ./Ouput/Rice/
```

**Dataset folder :**

It contains count tables with filename extension "fst" of Arabidopsis and Rice WT samples in this study.

First, all Degradome dataset were processed and normalized (TP40M) as described in this study.
Next, processed libraries were mapped to each transcript which concatenated with 100 nt of the corresponding transcript upstream and downstream sequences from the genome. 
Finally, the mapping results were summarized into the 5'P end count table (.tsv; tab-delimited file) and tsv files were coverted to fst format via R package "fst".

5'P end count table contains three columns of Reference, Position and Read (Dataset/example_count_table_cDNAU100D100.tsv).
<pre>
Frist  column (Reference, character) : transcript ID (eg. AT1G01010.1)

Second column (Position, number)     : 1-based position corresponds to the frist base of 100 nt of upstream seqeuence from the annotated transcript start site.

Third  column (Read, number)         : normalized abundance of 5'P end
</pre>

Execute R script "Tsv2Fst.R" will covert count table from tsv format to fst format, and it requires R package data.table and fst.
```
Rscript Tsv2Fst.R <Path_of_folder_contains_tsv_file>
```

**Reference folder :**

It contains annotation files
1. Annotation files are from TAIR(Arabidopsis;www.arabidopsis.org) and MSU7(Rice;rice.plantbiology.msu.edu)

2. tpc115485Supp-WT_vs_lba1upf3-1.tsv is from "(B)WT_vs_lba1upf3-1" sheet of Supplemental Data Set 1 in Drechsel et al. Plant Cell (2013).

3. Representative_gene_model_file (Representative_gene_model_With_annotation.tsv)

Representative_gene_model_file is merge files in Reference folder by R script "Representaive_Ann.R".

2019.11.30 written by Bo-Han Hou
