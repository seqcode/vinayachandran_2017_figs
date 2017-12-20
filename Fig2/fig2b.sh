set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files_b
TAR=GSE98573_RAW.tar

if [ ! -d $TAB_DIR ]
	then 
		mkdir $TAB_DIR	
		mv $TAR $TAB_DIR
		cd $TAB_DIR
		tar xvf $TAR
		rm $TAR
		gunzip *.gz
		cd ..
fi

GFF=ALL-RP-SAGA-TFIID-SUT-CUT-XUT_TSS-TES-MID_sortedby_geneLength.gff

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 -o $CDT_DIR $NORM_DIR ../shared_files/$GFF
fi

python ../scripts/sort_cdt_by_given_file.py -o 2 $CDT_DIR ../shared_files/$GFF 
