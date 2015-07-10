__author__ = 'krier'
from bein import *
from bein.util import *
import sys, os, re, json, shutil, gzip, tarfile, bz2, pickle, urllib, time
from bbcflib.common import set_file_descr

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

@program
def copy_here(file):
    shortname=os.path.basename(file)
    return {'arguments': ["cp"]+[file]+["."],
            'return_value': shortname}


@task
def copying(ex,files):
	for file in files.keys():
		print files[file]
		print file
		name="_".join(["Lib",str(file),"bwa.bam"])
		indexname=name+".bai"
		alignment=copy_here(ex,files[file][0])
		index=copy_here(ex,files[file][2])
		metric=copy_here(ex,files[file][3])
		metricname=name+".metrics"
		print name
		print indexname
        print metricname
		ex.add(alignment,alias=name,description=name)
		ex.add(index,alias=indexname,description=indexname,associate_to_filename=alignment,template="%s.bai")
		ex.add(metric,alias=metricname,description=metricname,associate_to_filename=alignment,template="%s.metrics")


libs=range(23,24)
iterate=str(1)

runs=libs

files={}

for run in runs:
	runname="_".join(["Lib",str(run),"bwa"])+".bam"
	files[run] = [f for f in listdir_fullpath('/data/leukemia_analysis/YmG3UHDpdo8pkhEb28sn/') if re.match(r'.*%s'%runname,f,re.IGNORECASE)]
	if "_R2_" in files[run][0]:
		files[run]=files[run][::-1]
	print files[run]

M=MiniLIMS("/data/leukemia_data/leukemia_debug")

copied=copying(M,files)