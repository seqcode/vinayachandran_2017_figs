set -e

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=(50519 50520 50521 50522 50523 50524 50525 50526 50527 50528 50529 50530)
TAB_DIR=tab_files_a

if [ ! -d $TAB_DIR ]
	then
		mkdir $TAB_DIR

		for ID in "${IDS[@]}"
		do
			cp $ALL_TAB/$ID"sacCer3".tab $TAB_DIR
		done
fi

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

OUT_DIRS=(RP SAGA TFIID CUTs SUTs XUTs)
GFFS=(RP_137_genes_TSS_Xu_2009.gff SAGA_TSS_Xu_2009_ORF_Ts_V64.gff TFIID_TSS_Xu_2009_ORF_Ts_V64.gff Xu_2009_CUTs_TSS_ONLY_V64.gff Xu_2009_SUTs_TSS_ONLY_V64.gff van_Dijk_2011_XUTs_V64_TSS_ONLY.gff)

for i in `seq 0 $((${#OUT_DIRS[@]}-1))`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}

	if [ ! -e ../shared_files/$GFF ]
		then
			wget 	#TODO: add url
			mv $GFF.gz ../shared_files
			gunzip ../shared_files/$GFF.gz
	fi

	if [ ! -d $OUT_DIR ]
		then
			python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR tab_files_a/Normalized_tab_files ../shared_files/$GFF
	fi
	python ../scripts/composite_plots.py -w 20 --shaded $OUT_DIR
done
