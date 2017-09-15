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
SORT_FILE=ALL-RP-SAGA-TFIID-SUT-CUT-XUT_TSS-TES-MID_sortedby_geneLength.gff

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 tab_files_b/Normalized_tab_files $SORT_FILE
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR $SORT_FILE
