set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 tab_files/Normalized_tab_files  Hsf1-union-Xu-TSS-divergent-upstream.gff
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR Hsf1-union-Xu-TSS-sortby-distance.gff
