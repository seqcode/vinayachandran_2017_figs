set -e

NORM_DIR=tab_files/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=CDT
SORT_FILE=../shared_files/ALL-RP-SAGA-TFIID-SUT-CUT-XUT_TSS-TES-MID_sortedby_geneLength.gff

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 -o $CDT_DIR $NORM_DIR $SORT_FILE
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR $SORT_FILE 
