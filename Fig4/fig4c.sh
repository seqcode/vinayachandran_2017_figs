set -e

sh ../scripts/get_chrom_sizes.sh

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=(50428 50429)
TAB_DIR=tab_files_c

if [ ! -d $TAB_DIR ]
	then
		mkdir $TAB_DIR

		for ID in "${IDS[@]}"
		do
			cp $ALL_TAB/$ID"sacCer3".tab $TAB_DIR
		done
fi

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

for TAB_FILE in $NORM_DIR/*.tab
do
	for ORIENTATION in upstream downstream
	do
		GFF="Hsf1-union-Xu-TSS-divergent-"$ORIENTATION.gff
		if [ ! -e $GFF ]
		then
			wget 	#TODO: add url
			gunzip $GFF.gz
		fi		

		python ../scripts/extract_tag_occupancy.py $TAB_FILE $GFF ../shared_files/sacCer3.chrom.sizes ${TAB_FILE%.*}$ORIENTATION.txt 50 50
	done
done

python fig4c.py
