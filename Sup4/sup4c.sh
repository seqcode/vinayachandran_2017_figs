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

python ../scripts/quantile_norm_singlebase_bin.py tab_files_c $CHROM_INFO

OUT_DIR=_CDT
if [ ! -d $OUT_DIR ]
	then
		mkdir $OUT_DIR
fi

python ../scripts/map_shifted_tags_to_ref.py -u 450 -d 450 -o $OUT_DIR tab_files_c/Normalized_tab_files $GFF	#TODO: where is gff
python ../scripts/composite_plots.py -w 20 --shaded $OUT_DIR
