set -e

NORM_DIR=tab_files_a/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_a ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=a_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR ../shared_files/TFIID_Plus_One_nucleosome_dyad.gff

fi

python ../scripts/composite_plots.py -w 20 --shaded $CDT_DIR
