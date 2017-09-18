set -e

WD=$PWD

SHARED_FILES=../shared_files
if [ ! -d $SHARED_FILES ]
	then
		mkdir $SHARED_FILES
fi

CHROM_INFO=../shared_files/sacCer3.chrom.sizes
if [ ! -e $CHROM_INFO ]
	then
		cd $SHARED_FILES
		wget http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
		cd $WD
fi

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py tab_files/Normalized_tab_files ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff
fi

python ../scripts/quantile_norm_singlebase_bin.py tab_files $CHROM_INFO --stranded
python ../scripts/Composite_for_many_factors_one_plot.py -w 20 tab_files/Normalized_tab_files 
