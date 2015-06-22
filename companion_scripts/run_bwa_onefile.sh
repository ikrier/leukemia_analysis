#!/bin/bash
<<"COMMENT"
The code should do the following :
Run alignment using bwa on the input fastq files for read1 and read2 in that order
then sort the file
then apply picard to mark duplicates in the file

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
COMMENT

usage ()
{
  echo 'Usage : run_bwa_onefile.sh <read1.fastq> <read2.fastq> <outputfile.srt.bam>'
  exit
}

if [ "$#" -ne 3 ]
then
  usage
fi

echo "Running bwa mem on" $1 "and" $2 "to create "$3

bwa mem -t 4 -M /data/genomes/Broadhs37/hs37d5.fa.gz $1 $2 2>$3.out|samtools view -bSu - |samtools sort - -f $3

samtools index $3

java -jar /data/software/picard/dist/picard.jar MarkDuplicates I=$3 O=$3temp M=$3.metrics &>>$3.out

mv $3temp $3

samtools index $3

cat $3.out
