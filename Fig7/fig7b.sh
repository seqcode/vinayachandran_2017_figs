set -e

SORT_FILE=sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt

if [ ! -e ../shared_files/$SORT_FILE ]
	then
		wget 	#TODO: add url
		mv $SORT_FILE.gz ../shared_files
		gunzip ../shared_files/$SORT_FILE.gz
fi

GFF=Xu_2009_RP-SAGA-TFIID_ONLY_TSS_ONLY_V64.gff

if [ ! -e ../shared_files/$GFF ]
	then
		wget 	#TODO: add url
		mv $GFF.gz ../shared_files
		gunzip ../shared_files/$GFF.gz
fi

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=(51829 51830 51831 56422 56423 56424 56428 56429 56430)
TAB_DIR=tab_files_b

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

CDT_DIR=b_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $CDT_DIR $NORM_DIR ../shared_files/$GFF
fi

python ../scripts/sort_cdt_by_given_file.py $CDT_DIR ../shared_files/$SORT_FILE
