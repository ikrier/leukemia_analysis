args <- commandArgs(trailingOnly = TRUE)

libname=args[1]
diamtablename=args[2]
fastqfile1=args[3]
fastqfile2=args[4]
bamname=args[5]
dupname=args[6]
regionsname=args[7]
ampliconsname=args[8]

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
