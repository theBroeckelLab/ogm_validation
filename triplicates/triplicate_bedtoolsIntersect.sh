#!/bin/bash

while getopts ":a:" opt; do
  case $opt in
    a)
      echo "-a was triggered, Parameter: $OPTARG" >&2
	name=${OPTARG}
	file1="${name}"
	file2="${name}_2"
	file3="${name}_3"
	# u means union file of the two listed
	u12="_12_u"
	# fr means fractional and recipricol parameters were set
	fr12="${name}_12_fr"
	u13="_13_u"
	fr12="_12_fr"
	#final files are the final union of 3 files , fr means f and r parameters were used
	finalu="_finalu"
	finalfr="_finalfr"
	# v means the unique calls to file
	v="_v"
	# vfr is unique calls to file with f and r parameters 
	vfr="_vfr"



	
	bedtools intersect -a ${file1}.bed -b ${file2}.bed ${file3}.bed -v  -header > ${file1}_v.bed
	bedtools intersect -a ${file2}.bed -b ${file1}.bed ${file3}.bed -v  -header > ${file2}_v.bed
	bedtools intersect -a ${file3}.bed -b ${file1}.bed ${file2}.bed -v  -header > ${file3}_v.bed
	bedtools intersect -a ${file1}.bed -b ${file2}.bed ${file3}.bed -v -f 0.9 -r -header > ${file1}_vfr.bed
	bedtools intersect -a ${file2}.bed -b ${file1}.bed ${file3}.bed -v -f 0.9 -r -header > ${file2}_vfr.bed
	bedtools intersect -a ${file3}.bed -b ${file1}.bed ${file2}.bed -v -f 0.9 -r -header > ${file3}_vfr.bed
	bedtools intersect -a ${file1}.bed -b ${file2}.bed -u  -header > ${file1}_12u.bed
	bedtools intersect -a ${file1}_12u.bed -b ${file3}.bed -u  -header > ${file1}_finalu.bed
	bedtools intersect -a ${file1}.bed -b ${file3}.bed -u  -header > ${file1}_13u.bed
	bedtools intersect -a ${file2}.bed -b ${file3}.bed -u  -header > ${file1}_23u.bed
	bedtools intersect -a ${file1}.bed -b ${file2}.bed -u -f 0.9 -r  -header > ${file1}_12fr.bed
	bedtools intersect -a ${file1}_12fr.bed -b ${file3}.bed -u -f 0.9 -r -header > ${file1}_finalfr.bed
	bedtools intersect -a ${file1}.bed -b ${file3}.bed -u -f 0.9 -r -header > ${file1}_13fr.bed
	bedtools intersect -a ${file2}.bed -b ${file3}.bed -u -f 0.9 -r  -header > ${file1}_23fr.bed

	bedtools intersect -a ${file1}_12fr.bed -b ${file3}.bed -v -f 0.9 -r -header > ${file1}_12vfr.bed
	bedtools intersect -a ${file1}_13fr.bed -b ${file2}.bed -v -f 0.9 -r -header > ${file1}_13vfr.bed
	bedtools intersect -a ${file1}_23fr.bed -b ${file1}.bed -v -f 0.9 -r  -header > ${file1}_23vfr.bed


	
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done
