#DianaOaxaca 231123
#QC filter, delete nextera adapter, Ns and low quality

mkdir -p results/02.trimgalore
out="results/02.trimgalore"
FASTQ=$(ls data/*.fastq.gz | sed 's/\(.*\)_.*/\1/' | sed 's/data\///' | sort -u)
declare -a FASTQs=("$FASTQ")

date

for FILE in ${FASTQs[@]}; do

	echo -e "Run trimgalore to $FILE sample"
	TRIMline='trim_galore --fastqc --length 250 -q 15 -j 12 --gzip  --paired data/'$FILE'_1.fastq.gz  data/'$FILE'_2.fastq.gz -o results/02.trimgalore/'$FILE'_trimgalore'
	echo -e $TRIMline "\n"
	$TRIMline
done

date
