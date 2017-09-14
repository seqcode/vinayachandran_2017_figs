set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files_a 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 tab_files_a/Normalized_tab_files sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt

fi

python ../scripts/composite_plots_shaded.py -w 21 $CDT_DIR
