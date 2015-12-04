#This is the script that should be called to map all the reads with bein.
""" The code should do the following :
bowtie2 [options]* -x <bowtie index prefix> {-1 <mate 1 file> -2 <mate 2 file>  [-S <sam output file>]
So we need the following parts :
- A function which submits the files for QC using fastqc or similar
- A function which takes a group of files in argument and submits in parallel the executions
- A function which provides a mapping quality report
- A function which provides a report on the coverage of the different regions we're examining,
taking as argument whether this is a design 2 or design 3 experiment depending.
Everything should be managed by bein and the working directory will be chosen as leukemia_data for the bein files

=======
License
=======
This code is released under the GNU General Public License 3.0. A copy
of this license is in the LICENSE.txt file.
copyright Irina Krier 2015
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=======
Attributions
=======
Parts of this code were adapted from bbcflib <https://github.com/bbcf/bbcflib> and such functions are indicated through comments
"""

from bein import *
from bein.util import *
import sys, os, re, json, shutil, gzip, tarfile, bz2, pickle, urllib, time
from bbcflib.common import set_file_descr


#Modified from bbcflib:
@program
def fastqc(fastqfile,outdir=".",options=None):
    outfile = re.sub(".fastq.gz","",os.path.basename(fastqfile))+'_fastqc.zip'
    if not(isinstance(options,list)): options = []
    if outdir and os.path.isdir(outdir):
        outfile = os.path.join(outdir,outfile)
        options += ["--outdir",outdir]
    return {'arguments': ["fastqc","--noextract"]+options+[fastqfile],'return_value': outfile}

@program
def bwa(fastqfiles,fileout):
    return {'arguments': ["/data/leukemia_analysis/companion_scripts/run_bwa_onefile.sh"]+fastqfiles+[fileout],
            'return_value': fileout}
@program
def bwa_notall(fastqfiles,fileout):
    return {'arguments': ["/data/leukemia_analysis/companion_scripts/run_bwa_onefile_notall.sh"]+fastqfiles+[fileout],
            'return_value': fileout}

@program
def markduplicates(fileout):
    name=fileout+".metrics"
    return {'arguments': ["/data/leukemia_analysis/companion_scripts/run_bwa_onefile_justmarkduplicates.sh"]+[fileout],
            'return_value': name}



@program
def cutadapt(fastqfiles,suffix,rm):
	name1=os.path.splitext(os.path.splitext(os.path.basename(fastqfiles[0]))[0])[0]+"."+suffix+".fastq.gz"
        name2=os.path.splitext(os.path.splitext(os.path.basename(fastqfiles[1]))[0])[0]+"."+suffix+".fastq.gz"
	return {'arguments': ["/data/leukemia_analysis/companion_scripts/cutadapt.sh"]+fastqfiles+[suffix]+[rm],'return_value': [name1 , name2]}

@program
def qcreport(libname,fastqfiles,bamfile):
	diamictable="/data/leukemia_analysis/diamtable.csv"
	if(float(libname)<48):
		targets="/data/generate_QC_reports/2ndDesignRegions.bed"
		amplicons="/data/generate_QC_reports/15189-1368439012_Amplicons.bed"
	else:
		targets="/data/generate_QC_reports/3rdDesignRegions.bed"
		amplicons="/data/generate_QC_reports/15189-1405499558_Amplicons.bed"
	libname=str(libname)
	outfile="QCreport_"+libname+".pdf"
	return {'arguments': ["Rscript","/data/leukemia_analysis/QC_report/run_sweave_libname.R"]+[libname]+[diamictable]+fastqfiles+bamfile+[targets]+[amplicons], 'return_value': outfile}

@program
def mergefiles(bamfiles,outname):
	return {'arguments': ["samtools","merge"]+[outname]+bamfiles,'return_value': outname}

def add_file_fastqc(execution,filename, description="", alias="None"):
    execution.add(filename,description=description,alias=alias)

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

#I should use an association between the two files, with a pattern using _R1.fastq.gz and _R2.fastq.gz template
@task
def trim_adapt(ex,files,suffix):
	for file in files.keys():
		print files[file]
		print file
		name1="Lib_"+str(file)+"_R1_"+suffix+".fastq.gz"
		name2="Lib_"+str(file)+"_R2_"+suffix+".fastq.gz"
		outputs=cutadapt(ex,files[file],suffix,"5")
		ex.add(outputs[0],alias=name1,description=name1)
		ex.add(outputs[1],alias=name2,description=name2)
	return {"test":"test"}
#To use a file you can use execution.use()

@task
def align_bwa(ex,files):
	for file in files.keys():
		print files[file]
		print file
		name="_".join(["Lib",str(file),"bwa.bam"])
		indexname=name+".bai"
		alignment=bwa(ex,files[file],name)
		index=alignment+".bai"
		metric=alignment+".metrics"
		metricname=name+".metrics"
		print name
		print indexname
		ex.add(alignment,alias=name,description=name)
		ex.add(index,alias=indexname,description=indexname,associate_to_filename=alignment,template="%s.bai")
		ex.add(metric,alias=metricname,description=metricname,associate_to_filename=alignment,template="%s.metrics")

@task
def align_bwa_notall(ex,files):
        for file in files.keys():
                print files[file]
                print file
                name="_".join(["Lib",str(file),"bwa.bam"])
                indexname=name+".bai"
		print "flag1"
                alignment=bwa_notall(ex,files[file],name)
                index=alignment+".bai"
		print "flag2"
                metric=markduplicates(ex,alignment)
		print "flag3"
                metricname=name+".metrics"
                print name
                print indexname
		print metric
		print metricname
                ex.add(alignment,alias=name,description=name)
                ex.add(index,alias=indexname,description=indexname,associate_to_filename=alignment,template="%s.bai")
                ex.add(metric,alias=metricname,description=metricname,associate_to_filename=alignment,template="%s.metrics")


@task
def run_qc_report(ex,fastqfiles,bamfiles):
	for file in fastqfiles.keys():
		print fastqfiles[file]
		print bamfiles[file]
		report=qcreport(ex,file,fastqfiles[file],bamfiles[file])
		reportname="QCreport_Lib"+str(file)+".pdf"
		ex.add(reportname,alias=reportname,description="QC report for "+str(file))

#Add new tasks running the variant callers, i.e. varscan, mutect and GATK
