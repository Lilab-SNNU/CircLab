#!/bin/bash
circ_bwa=$1
circ_star=$2
ciri=$3
find_circ=$4
output_bed1=$5
output_bed2=$6
output_bed3=$7
output_bed4=$8

# 先清空目标文件以确保开始时它是空的
> "$output_bed1"
> "$output_bed2"
> "$output_bed3"
> "$output_bed4"

cut -f 1,2,3,6 "$circ_bwa" > "${circ_bwa}_name0.txt"
cut -f 1,2,3,6 "$circ_star" > "${circ_star}_name0.txt"
cut -f 2,3,4,11 "$ciri" > "${ciri}_name0.txt"

tail -n +2 "${ciri}_name0.txt" > "${ciri}_name1.txt"
awk '{print $1, $2, $3,"",".", $4, $5}' OFS="\t" "${circ_bwa}_name0.txt" > "${circ_bwa}_name2.txt"
awk '{print $1, $2, $3,"",".", $4, $5}' OFS="\t" "${circ_star}_name0.txt" > "${circ_star}_name2.txt"
awk '{print $1, $2, $3,"",".", $4, $5}' OFS="\t" "${ciri}_name1.txt" > "${ciri}_name2.txt"

# 在${circ_bwa}_name.txt中添加一个新的字段，其内容是id，并直接输出到目标文件
awk -v id="$sample" 'BEGIN{FS=OFS="\t"} { if ($6 == "+") {$4 = "forward\t" $4} else if ($6 == "-") {$4 = "reverse\t" $4} print $0, id }' "${circ_bwa}_name2.txt" > "${circ_bwa}_name.txt"
awk -v id="$sample" 'BEGIN{FS=OFS="\t"} { if ($6 == "+") {$4 = "forward\t" $4} else if ($6 == "-") {$4 = "reverse\t" $4} print $0, id }' "${circ_star}_name2.txt" > "${circ_star}_name.txt"
awk -v id="$sample" 'BEGIN{FS=OFS="\t"} { if ($6 == "+") {$4 = "forward\t" $4} else if ($6 == "-") {$4 = "reverse\t" $4} print $0, id }' "${ciri}_name2.txt" > "${ciri}_name.txt"

cat "${circ_bwa}_name.txt" >> "$output_bed1"
cat "${circ_star}_name.txt" >> "$output_bed2"
cat "${ciri}_name.txt" >> "$output_bed3"
cat "${find_circ}" >> "$output_bed4"

rm "${circ_bwa}_name.txt" "${circ_bwa}_name0.txt" "${circ_bwa}_name2.txt"
rm "${circ_star}_name.txt" "${circ_star}_name0.txt" "${circ_star}_name2.txt"
rm "${ciri}_name.txt" "${ciri}_name0.txt" "${ciri}_name1.txt" "${ciri}_name2.txt"
