import vcf
import os
import subprocess

vcf_reader=vcf.Reader(open("all_merged.vcf.gz","r"))
		
vcf_writer = vcf.Writer(open("all_merged.flt.vcf", 'w'), vcf_reader)
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
	if float(db.get("g1000",0))<0.05 and float(db.get("cg69",0))<0.05 and not record.FILTER and maxAF>0.1:
		vcf_writer.write_record(record)
	
vcf.Writer.close(vcf_writer)

