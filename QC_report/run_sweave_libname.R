setwd("/data/generate_QC_reports/temp_temp_dir/")
costable=read.csv("cosmtable.csv",header=F)

fastqdir="/data/fastq_gwendal/"
libname="49"
cosmname=costable$V2[as.numeric(libname)]
if(cosmname=="NONE"){cosmname="Aucun"}

regionsname="../3rdDesignRegions.bed"
ampliconsname="../15189-1405499558_Amplicons.bed"
#For the 2nd design : 
# regionsname="../2ndDesignRegions.bed"
# ampliconsname="../15189-1368439012_Amplicons.bed"

reportname=paste("Report",libname,"_draft",sep="")

Sweave(paste(reportname,".Rnw",sep=""))
tools::texi2pdf(paste(reportname,".tex",sep=""))
