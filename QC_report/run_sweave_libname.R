# Here are the libraries to load for this to work :
library(xtable)
library("qrqc")
library("RColorBrewer")
library("RColorBrewer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("rtracklayer")
library(RColorBrewer)

#In addition the program requires the following programs to be installed :
# Picard
# Samtools
# Bedtools

#This program can be run providing a number of arguments read from the command line :

args <- commandArgs(trailingOnly = TRUE)

libname=args[1] #Library name
diamtablename=args[2] #A correspondance table between library name and DIAMIC sample number
fastqfile1=args[3] #The first fastq file for paired sequencing
fastqfile2=args[4] #The second
bamname=args[5] #The bam library
dupname=args[6] #The bam library with duplicates marked
regionsname=args[7] #The bed file of the target regions
ampliconsname=args[8] #The bed file of the amplicons

# Example command to run the script :
# Rscript run_sweave_libname.R 1 table_samples.csv Lib_1_R1.fastq Lib_1_R2.fastq Lib_1.bam Lib_1.dups.bam target_regions.bed amplicons_regions.bed
# "diamtable" would be comma-separated values with the first being the library name and the second the sample ID.
# Working example :
# Rscript /data/leukemia_analysis/QC_report/run_sweave_libname.R 49 /data/leukemia_analysis/diamtable.csv /data/TruSight_analysis/BM-cyto-1_S1_R1.fastq.gz  /data/TruSight_analysis/BM-cyto-1_S1_R2.fastq.gz /data/TruSight_analysis/Lib_49_realigned_cutends.srt.bam /data/TruSight_analysis/Lib_49_realigned_cutends.srt.dups.bam /data/TruSight_analysis/Targets_disjoint_noprimers.trusight.bed /data/TruSight_analysis/trusight-myeloid-amplicon-track.bed 

for(i in 1:7)
{
  cat(args[i],"\n")
}

.libPaths( c( .libPaths(), "/home/krier/R/x86_64-pc-linux-gnu-library/3.1/") )

system(command = paste("mkdir qcreportdir",libname,sep="_"))
setwd(paste("qcreportdir",libname,sep="_"))

diamtable=read.csv(diamtablename,header=F)

diamname=diamtable$V2[as.numeric(libname)]
if(diamname=="NONE"){diamname="Aucun"}

#For the 3rd design :
#regionsname="../3rdDesignRegions.bed"
#ampliconsname="../15189-1405499558_Amplicons.bed"
#For the 2nd design : 
# regionsname="../2ndDesignRegions.bed"
# ampliconsname="../15189-1368439012_Amplicons.bed"

Sweave("/data/leukemia_analysis/QC_report/QCreport.Rnw")
tools::texi2pdf("QCreport.tex")

system(command=paste("mv QCreport.pdf ../QCreport_Lib",libname,".pdf",sep=""))
