set -e

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

SORT_FILE=ALL-RP-SAGA-TFIID-SUT-CUT-XUT_TSS-TES-MID_sortedby_geneLength.gff

if [ ! -e ../shared_files/$SORT_FILE ]
	then
		wget 	#TODO: add url
		mv $SORT_FILE.gz ../shared_files
		gunzip ../shared_files/$SORT_FILE.gz
fi

IDS=(50416 50513 50501 50525 53319 50428 50531 50519 54513 51825 53803 50417 50514 50502 50526 53320 50429 50532 50520 54514 53407 53804)
TAB_DIR=tab_files

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

CDT_DIR=CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 -o $CDT_DIR $NORM_DIR ../shared_files/$SORT_FILE
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR ../shared_files/$SORT_FILE 
