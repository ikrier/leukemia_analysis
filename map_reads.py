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

M = MiniLIMS("/data/leukemia_data/minilims")

with execution(M) as ex:
   touch(ex, "boris")
   print ex.working_directory

@program
def fastqc(fastqfile,outdir=None,options=None):
    """Binds ``fastqc`` (`<http://www.bioinformatics.bbsrc.ac.uk/>`_) which generates a QC report of short reads present in the fastq file.
    """
    outfile = re.sub(".fastq","",os.path.basename(fastqfile))+'_fastqc.zip'
    if not(isinstance(options,list)): options = []
    if outdir and os.path.isdir(outdir):
        outfile = os.path.join(outdir,outfile)
        options += ["--outdir",outdir]
    return {'arguments': ["fastqc","--noextract"]+options+[fastqfile],'return_value': outfile}

def run_fastqc( ex, job):
    """
    Returns the name of the report file.
    """
    futures = {}
    descr = {'step':'qc','groupId':0,'type':'zip'}
    for gid,group in job.groups.iteritems():
        futures[gid] = {}
        for rid,run in group['runs'].iteritems():
            if isinstance(run,tuple):
                futures[gid][rid] = (fastqc.nonblocking(ex,run[0],via=via),
                                     fastqc.nonblocking(ex,run[1],via=via))
            else:
                futures[gid][rid] = fastqc.nonblocking(ex,run,via=via)
    for gid,group in job.groups.iteritems():
        descr['groupId'] = gid
        for rid,run in group['runs'].iteritems():
            rname = group['name']
            if len(group['runs'])>1:
                rname += "_"
                rname += group['run_names'].get(rid,str(rid))
            if isinstance(run,tuple):
                qcreport = futures[gid][rid][0].wait()
                if os.path.exists(qcreport):
                    ex.add( qcreport,
                            description=set_file_descr(rname+"_R1_fastqc.zip",**descr) )
                qcreport = futures[gid][rid][1].wait()
                if os.path.exists(qcreport):
                    ex.add( qcreport,
                            description=set_file_descr(rname+"_R2_fastqc.zip",**descr) )
            else:
                qcreport = futures[gid][rid].wait()
                if os.path.exists(qcreport):
                    ex.add( qcreport,
                            description=set_file_descr(rname+"_fastqc.zip",**descr) )
    return None