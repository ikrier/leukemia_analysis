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
from map_reads import *

libs=range(int(sys.argv[1]),int(sys.argv[2]+1))
iterate=sys.argv[3]

runs=libs

files={}

for run in runs:
	runname="_".join(["Lib",str(run)])
	files[run] = [f for f in listdir_fullpath('/data/fastq_gwendal/') if re.match(r'.*%s_'%runname,f,re.IGNORECASE)]
	if "_R2_" in files[run][0]:
		files[run]=files[run][::-1]
	print files[run]

M=MiniLIMS("/data/leukemia_data/leukemia_pipeline")

suffix="adaptedtrimmed5"
trimming=trim_adapt(M,files,suffix)
with execution(M) as ex:
	add_pickle(ex,trimming,description="object for trimmed files info",alias="trimming"+iterate)
#to use it after again : use_pickle(M,"trimming1")

trimmedfiles={}

for run in runs:
	trimname1="Lib_"+str(run)+"_R1_"+suffix+".fastq.gz"
        trimname2="Lib_"+str(run)+"_R2_"+suffix+".fastq.gz"
	trimmedfiles[run]=[M.path_to_file(trimming["files"][trimname1]), M.path_to_file(trimming["files"][trimname2])]

aligning=align_bwa(M,trimmedfiles)
with execution(M) as ex:
	add_pickle(ex,aligning,description="object for aligned files info",alias="aligned"+iterate)

bamfiles={}

for run in runs:
	bamfiles[run]=[M.path_to_file(aligning["files"]["Lib_"+str(run)+"_bwa.bam"])]

qcreporting=run_qc_report(M,trimmedfiles,bamfiles)
with execution(M) as ex:
	add_pickle(ex,qcreporting,description="object for qc reporting files info",alias="qcreporting"+iterate)
