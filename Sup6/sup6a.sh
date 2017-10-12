set -e

NORM_DIR=tab_files_a/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_a ../shared_files/sacCer3.chrom.sizes
fi

OUT_DIR=a_CDT
if [ ! -d $OUT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 30 -d 100 -o $OUT_DIR $NORM_DIR ../shared_files/TATA_consensus_anti_Max_sort.bed
fi

python ../scripts/composite_plots.py -w 20 $OUT_DIR
