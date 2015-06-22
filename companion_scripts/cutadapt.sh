#!/bin/bash
<<"COMMENT"
 The code should do the following :
Remove Illumina adaptors from the fastq files.
Remove the first 5bp of each read as they contain the restriction enzyme recognition sequences and therefore no mismatch will be mapped there.

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
  echo 'Usage : ./cutadapt.sh <fastq1> <fastq2> <suffix> removestart'
  exit
}

if [ "$#" -ne 4 ]
then
  usage
fi

file1=$(basename $1)
file2=$(basename $2)

name1=${file1%%.*}
name2=${file2%%.*}

echo $name1
echo $name2

echo "Running cutadapt on" $1 "and" $2 "to create "$name1.$3.fastq.gz and $name2.$3.fastq.gz


cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $1 $2 -o $name1.$3.fastq.gz -p $name2.$3.fastq.gz -m 20 &>$3.out

#for optional trimming :
cutadapt -u $4 -o $name1.$3.trim5.fastq.gz $name1.$3.fastq.gz &>>$3.out
cutadapt -u $4 -o $name2.$3.trim5.fastq.gz $name2.$3.fastq.gz &>>$3.out

mv $name1.$3.trim5.fastq.gz $name1.$3.fastq.gz
mv $name2.$3.trim5.fastq.gz $name2.$3.fastq.gz

cat $3.out
