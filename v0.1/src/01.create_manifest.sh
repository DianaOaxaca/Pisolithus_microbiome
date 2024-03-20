#!/usr/bin/bash
#DianaOaxaca 200224
#Create manifest file

cd data/
id=$(ls *.gz | sed 's/_.*//g' | sort -u)

echo -e "sample-id"",""absolute-filepath"",""direction" > manifest.csv

for sample in $id; do
	source=$(awk -F"\t" -v sample="$sample" '$2 == sample {print $1}' metadata.tsv)
        r1=$(echo -e $sample"_1.fastq.gz"",""forward")
        r2=$(echo -e $sample"_2.fastq.gz"",""reverse")
        path=$(pwd)
        echo -e $source","$path"/"$r1 >> manifest.csv
        echo -e $source","$path"/"$r2 >> manifest.csv
done
