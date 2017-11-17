set -e

sh ../scripts/get_chrom_sizes.sh

GFF=TFIID-nochange_TSS_Xu_2009_ORF_PURE_SET.gff

if [ ! -e $GFF ]
	then
		wget 	#TODO: add url
		gunzip $GFF.gz
fi

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=()		#TODO: add IDs
TAB_DIR=tab_files

if [ ! -d $TAB_DIR ]
	then
		mkdir $TAB_DIR

		for ID in "${IDS[@]}"
		do
			cp $ALL_TAB/$ID"sacCer3".tab $TAB_DIR
		done
fi

#get TFIIH occupancy
python ../scripts/extract_tag_occupancy.py 19325sacCer3.tab $GFF ../shared_files/sacCer3.chrom.sizes TFIIH_occupancy.txt 100 100

#sort occupancy file
cat $GFF | awk '{print $9}' > genes.txt
paste genes.txt TFIIH_occupancy.txt > TFIIH_occupancy_by_gene.tsv
sort -k 2 -n TFIIH_occupancy_by_gene.tsv | awk '{print $1}' > sorted_TFIIH_genes.tsv
QUARTER_LENGTH=$(($(cat sorted_TFIIH_genes.tsv | wc -l)/4))
cat sorted_TFIIH_genes.tsv | head -$QUARTER_LENGTH > bottom_25percent.txt
cat sorted_TFIIH_genes.tsv | tail -$QUARTER_LENGTH > top_25percent.txt

#convert gene lists to GFF
python genes_from_gff.py $GFF bottom_25percent.txt
python genes_from_gff.py $GFF top_25percent.txt

#map to each GFF
python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o bottom25 tab_files bottom_25percent.gff
python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o top25 tab_files top_25percent.gff

#make composite plots
python ../scripts/composite_plots.py --normalize bottom25
python ../scripts/composite_plots.py --normalize top25
