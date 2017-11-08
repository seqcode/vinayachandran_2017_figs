set -e

TAB_DIR=tab_files_a

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=a_CDT
if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 100 -d 100 -o $CDT_DIR $NORM_DIR ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff
fi

python sup1a.py $CDT_DIR
