set -e

POLYA_GFF=polyA_sacCer3.gff
if [ ! -e ../shared_files/$POLYA_GFF ]
	then
		wget 	#TODO: add url
		mv $POLYA_GFF.gz ../shared_files
		gunzip ../shared_files/$POLYA_GFF.gz
fi

POLYT_GFF=polyT_sacCer3.gff
if [ ! -e ../shared_files/$POLYT_GFF ]
	then
		wget 	#TODO: add url
		mv $POLYT_GFF.gz ../shared_files
		gunzip ../shared_files/$POLYT_GFF.gz
fi

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=50502
TAB_DIR=tab_files_c

if [ ! -d $TAB_DIR ]
	then
		mkdir $TAB_DIR

		cp $ALL_TAB/$ID"sacCer3".tab $TAB_DIR
fi

POLYA_CDT_DIR=polyA_CDT

if [ ! -d $POLYA_CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $POLYA_CDT_DIR $TAB_DIR ../shared_files/$GFF

fi

python ../scripts/composite_plots.py $POLYA_CDT_DIR

POLYT_CDT_DIR=polyT_CDT

if [ ! -d $POLYT_CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $POLYT_CDT_DIR $TAB_DIR ../shared_files/$GFF

fi

python ../scripts/composite_plots.py $POLYT_CDT_DIR
