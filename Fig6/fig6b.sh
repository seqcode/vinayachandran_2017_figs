set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files_b 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 tab_files_b/Normalized_tab_files HS_activated_nuc1.gff
fi

python ../scripts/composite_plots_shaded.py -w 20 $CDT_DIR
