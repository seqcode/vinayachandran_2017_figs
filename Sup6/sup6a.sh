set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files_a 

OUT_DIR=_CDT
if [ ! -d $OUT_DIR ]
	then
		mkdir $OUT_DIR
fi

python ../scripts/map_shifted_tags_to_ref.py -u 30 -d 100 tab_files_a/Normalized_tab_files TATA_consensus_anti_Max_sort.bed
python ../scripts/Composite_plots_Vinesh\(1\).py -w 20 $OUT_DIR
