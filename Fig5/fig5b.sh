set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 tab_files/Normalized_tab_files  Yeast_plus_one_sacCer3.gff

fi

python ../scripts/composite_plots_shaded.py -w 20 $CDT_DIR
