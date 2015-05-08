import vcf
import os
import subprocess

""" The code should do the following :
Filter stuff
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

"""


filein="all_merged.vcf.gz"
fileout="all_merged.flt.vcf"

vcf_reader=vcf.Reader(open(filein,"r"))
		
vcf_writer = vcf.Writer(open(fileout, 'w'), vcf_reader)
for record in vcf_reader:
#We only take mutations which are present in less than 5% of the 1000 genomes and the cg69 genomes, that haven't been filtered by Sophia as being off target or anything else, and where the minor allele percentage is at least 10% because lower values may be irreproducible :
	DBinfo=record.INFO.get("DBXREF","")
	x=[]
	y=[]
	for item in DBinfo:
		if len(item.split(":"))>1:
			x.append(item.split(":")[0])
			y.append(item.split(":")[1])
		db=dict(zip(x,y))
	maxAF=0.1
	for sample in record.samples:
		if sample.data.AD:
			AF=map(float,sample.data.AD)[1]/sum(map(float,sample.data.AD))
			if AF>maxAF:
				maxAF=AF
	if float(db.get("g1000",0))<0.05 and float(db.get("cg69",0))<0.05 and maxAF>0.1 and (not record.FILTER or record.FILTER==["off_target"]):
		vcf_writer.write_record(record)
	
vcf.Writer.close(vcf_writer)

