set -e

NORM_DIR=tab_files_b/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_b ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT
SORT_FILE=../shared_files/sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $CDT_DIR $NORM_DIR $SORT_FILE
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR $SORT_FILE
