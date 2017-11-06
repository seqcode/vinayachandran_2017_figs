set -e

NORM_DIR=tab_files_b/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_b ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT
SORT_FILE=../shared_files/Xu_2009_ORF_TSS_TES_V64.gff	#TEST

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 -o $CDT_DIR tab_files_b/Normalized_tab_files $SORT_FILE
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR $SORT_FILE 
