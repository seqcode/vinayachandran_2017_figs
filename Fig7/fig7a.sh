set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files_a
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

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=a_CDT

#TODO: split by condition
if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 -o $CDT_DIR $NORM_DIR ../shared_files/RP_137_genes_TSS_Xu_2009.gff

fi

python ../scripts/composite_plots.py -w 21 --normalize $CDT_DIR
