set -e

NORM_DIR=tab_files/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR ../shared_files/sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt
