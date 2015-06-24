# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

plot_gene_profile=function(range,bamname,rmdupfile,targets,amplicons,title)
{
library("Gviz")
library("biomaRt")
library("BiocGenerics")
library("RColorBrewer")
ensembl67 <- useMart(host='may2012.archive.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL') #uses ensemble67 build ("may2012.archive.ensembl.org")
ensembl67 <- useDataset("hsapiens_gene_ensembl", mart=ensembl67)

pad=100
start(range)=start(range)-pad
end(range)=end(range)+pad

biomTrack <- BiomartGeneRegionTrack(genome = "hg19", biomart = ensembl67,
                                    chromosome = seqnames(range), stacking="dense",
                                    showId = TRUE, name= "Gene",start = start(range), end = end(range))
#if ensembl is not available for some reason :
#load(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txDB= TxDb.Hsapiens.UCSC.hg19.knownGene
#genestrack=GeneRegionTrack(txDB,showId=TRUE)
#Won't show the proper gene names but will work
cosmic=read.table("/data/leukemia_analysis/QC_report/COSMIC_Haematopoietic_and_lymphoid_tissue.gz",skip=1,header=T,sep="\t")
cosmic=GRanges(seqnames = cosmic$chromosome,ranges = IRanges(start=cosmic$grch37_start,end=cosmic$grch37_stop),
               strand="*",name=cosmic$mut_syntax_aa)
options(ucscChromosomeNames=FALSE)
cosmicTrack=AnnotationTrack(cosmic,name = "COSMIC Haematopoietic",stacking = "dense")
library(BSgenome.Hsapiens.UCSC.hg19)
options(ucscChromosomeNames=FALSE)
lengths=seqlengths(Hsapiens)
names(lengths)=sapply(names(lengths),function(x){substr(x,start = 4,stop = nchar(x))})
source("/data/leukemia_analysis/QC_report/strandedBamImport.R")
granges=strandedBamImport(bamname,range)
granges2=granges
granges2$both=NULL
rm(granges)
grangesu=strandedBamImport(rmdupfile,range)
granges2u=grangesu
granges2u$both=NULL
rm(grangesu)
gtrack <- GenomeAxisTrack()
dTrackGrange2=DataTrack(granges2,ylim=c(-17500,17500),name="Stranded coverage")
dTrackGrange2u=DataTrack(granges2u,ylim=c(-50,50),name="Unique coverage")
dTrackns <- DataTrack(range=bamname, genome="hg19", name="Total coverage",
                      window=-1, chromosome=seqnames(range),
                      stream=TRUE,ylim=c(0,35000))
displayPars(dTrackGrange2)$groups=c("+","-")
displayPars(dTrackGrange2)$col= c("#0099FF","#FF33CC")
displayPars(dTrackGrange2u)$groups=c("+","-")
displayPars(dTrackGrange2u)$col= c("#0099FF","#FF33CC")
targetsTrack=AnnotationTrack(targets,name="Targets")
displayPars(targetsTrack)$fill= "grey"
displayPars(targetsTrack)$col= "grey"
ampliconsTrack=AnnotationTrack(amplicons,stacking = "squish",name="Ampli- cons")
strand(ampliconsTrack)="*"
plotTracks(list(gtrack,dTrackns,dTrackGrange2,dTrackGrange2u,cosmicTrack,biomTrack,targetsTrack,ampliconsTrack),
           from =start(range), to = end(range),chromosome=seqnames(range),
           type="hist", col.histogram=NA,windowSize=0,baseline=0,sizes = c(1,1,1,1,1,0.5,0.5,0.5),
           main=title)
}