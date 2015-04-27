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

hg19_path="/data/genomes/hg19.bwa/hg19.fa"

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
def _bwa(fastqfiles,options=None):
    if "_R2_" in fastqfiles[0]:
        fastqfiles=fastqfiles[::-1]
    return {'arguments': ["bwa","mem", "-t 4"]+options+[hg19_path]+fastqfiles,
            'return_value': None}

def bwa(ex, fastqfiles, options=None,filename=None):
    if filename == None:
        filename = unique_filename_in()
    _bwa(ex, fastqfiles, options, stdout=filename)
    return filename



def add_file_fastqc(execution,filename, description="", alias="None"):
    execution.add(filename,description=description,alias=alias)

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

M = MiniLIMS("/data/leukemia_data/fastqreports")
#fileid=M.import_file(fastqfiles[0])

runs=range(1,85)

files={}

# for run in runs:
#     runname="_".join(["Lib",str(run)])
#     files[run] = [f for f in listdir_fullpath('/data/fastq_gwendal') if re.match(r'.*%s_'%runname,f,re.IGNORECASE)]
#     print files[run]

#files=dict((k, files[k]) for k in (1,2,3))

# with execution(M) as ex:
#     outdir="."
#     for file in files.keys():
#         print files[file]
#         # report=fastqc(ex,files[file][0],outdir=outdir,options=None)
#         # add_file_fastqc(ex,report,alias="_".join(["Lib",str(files.keys()[file-1]),"report_1"]))
#         # report=fastqc(ex,files[file][1],outdir=outdir,options=None)
#         # add_file_fastqc(ex,report,alias="_".join(["Lib",str(files.keys()[file-1]),"report_2"]))
#         futures=[fastqc.nonblocking(ex,f,outdir=outdir,options=None) for f in files[file]]
#         reports=[f.wait() for f in futures]
#         for report in reports:
#             i=str(report)[-17:-15]
#             add_file_fastqc(ex,report,alias="_".join(["Lib",str(files.keys()[file-1]),"report",i]))
#

with execution(M) as ex:
    for file in files.keys():
        print files[file]
        alignment=bwa(ex,files[file])
        name="_".join(["Lib",str(files.keys()[file-1]),"bwa.bam"])
        ex.add(alignment,alias=)