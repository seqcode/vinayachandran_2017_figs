set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files
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

SORT_FILE=sorted_list-RPG_SAGA-act-rep-no_TFIID-act-rep-no.txt

GFF=Xu_2009_RP-SAGA-TFIID_ONLY_TSS_ONLY_V64.gff

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR ../shared_files/$GFF
fi

python ../scripts/sort_cdt_by_given_file.py $CDT_DIR ../shared_files/$SORT_FILE
