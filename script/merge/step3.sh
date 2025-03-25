#!/bin/bash
ref=$1
filter1=$2
filter2=$3
filter3=$4
filter4=$5
circseq1=$6
circseq2=$7
circseq3=$8
circseq4=$9

#bedtools getfasta -fi ${ref} -bed "${filter1}" -fo "${circseq1}"0 -s
#bedtools getfasta -fi ${ref} -bed "${filter2}" -fo "${circseq2}"0 -s
#bedtools getfasta -fi ${ref} -bed "${filter3}" -fo "${circseq3}"0 -s
#bedtools getfasta -fi ${ref} -bed "${filter4}" -fo "${circseq4}"0 -s
bedtools getfasta -fi ${ref} -bed "${filter1}" -fo "${circseq1}" -s
bedtools getfasta -fi ${ref} -bed "${filter2}" -fo "${circseq2}" -s
bedtools getfasta -fi ${ref} -bed "${filter3}" -fo "${circseq3}" -s
bedtools getfasta -fi ${ref} -bed "${filter4}" -fo "${circseq4}" -s
#cd-hit-est -i "${circseq1}"0 -o "${circseq1}" -n 10
#cd-hit-est -i "${circseq2}"0 -o "${circseq2}" -n 10
#cd-hit-est -i "${circseq3}"0 -o "${circseq3}" -n 10
#cd-hit-est -i "${circseq4}"0 -o "${circseq4}" -n 10

#rm "${circseq1}"0 "${circseq2}"0 "${circseq3}"0 "${circseq4}"0
