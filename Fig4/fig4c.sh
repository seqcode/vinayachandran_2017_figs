set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files_c 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 tab_files_c/Normalized_tab_files  Hsf1-union-Xu-TSS-divergent-downstream.gff
fi

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 tab_files_c/Normalized_tab_files  Hsf1-union-Xu-TSS-divergent-upstream.gff
fi

python ../scripts/extract_tag_occupancy.py $FN 50 50	#TODO: what is fn??
