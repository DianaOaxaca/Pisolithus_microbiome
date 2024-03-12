#DianaOaxaca 291123
#Remove primers

mkdir -p results/03.cutadapt
out="results/03.cutadapt"
FASTQ=$(ls data/trim/*.gz | sed 's/_._val_.*.*//g' | sed 's/data\/trim\///g'  | sort -u)
declare -a FASTQs=("$FASTQ")

date

for FILE in ${FASTQs[@]}; do
	echo -e "Run cutadapt to $FILE sample"
        CUTADAPT='cutadapt --pair-filter any --no-indels -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC --discard-untrimmed -Z -j 40 -o '$out'/'$FILE'_1.fastq.gz -p '$out'/'$FILE'_2.fastq.gz data/trim/'$FILE'_1_val_1.fq.gz  data/trim/'$FILE'_2_val_2.fq.gz'
        echo -e $CUTADAPT  "\n"
        $CUTADAPT
done

date
