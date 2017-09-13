set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

TAB_DIR=tab_files
if [ ! -d $DIVERGENT_DIR ]
	then
		mkdir $DIVERGENT_DIR
fi

if [ ! -d $NONDIVERGENT_DIR ]
	then
		mkdir $NONDIVERGENT_DIR
fi

python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $DIVERGENT_DIR $TAB_DIR Hsf1-union-Xu-TSS-divergent-upstream.gff
python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $NONDIVERGENT_DIR $TAB_DIR Hsf1-union-Xu-TSS-sortby-distance.gff

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR Hsf1-union-Xu-TSS-sortby-distance.gff
