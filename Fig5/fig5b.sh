set -e

WD=$PWD

SHARED_FILES=../shared_files
if [ ! -d $SHARED_FILES ]
	then
		mkdir $SHARED_FILES
fi

CHROM_INFO=../shared_files/hg19.chrom.sizes
if [ ! -e $CHROM_INFO ]
	then
		cd $SHARED_FILES
		wget https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
		cd $WD
fi

python ../scripts/quantile_norm_singlebase_bin.py tab_files_b $CHROM_INFO

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 tab_files_b/Normalized_tab_files  Yeast_plus_one_sacCer3.gff

fi

python ../scripts/composite_plots.py -w 20 --shaded $CDT_DIR
