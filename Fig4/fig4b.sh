set -e

NORM_DIR=tab_files_b/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_b ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 400 -d 400 -o $CDT_DIR $NORM_DIR ../shared_files/Hsf1-union-midpoint.gff
fi

python ../scripts/composite_plots.py -w 20 --shaded --normalize $CDT_DIR
