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

FACTORS=(Ssl2 Spt3)
ALL_IDS=(50519,50520,50521,50522,50523,50524 50525,50526,50527,50528,50529,50530)

for i in `seq 0 $((${#FACTORS[@]}-1))`
do
	FACTOR=${FACTORS[$i]}
	IDS=${ALL_IDS[$i]}
	IDS_ARRAY=(${IDS//,/ })
	FACTOR_DIR=$TAB_DIR/$FACTOR

	if [ ! -d $FACTOR_DIR ]
		then
			mkdir $FACTOR_DIR

			for ID in "${IDS_ARRAY[@]}"
			do
				mv $TAB_DIR/*$ID"sacCer3".tab $FACTOR_DIR
			done
	fi

	NORM_DIR=$FACTOR_DIR/Normalized_tab_files

	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py $FACTOR_DIR ../shared_files/sacCer3.chrom.sizes
	fi

	GENE_CLASSES=(RP SAGA TFIID CUTs SUTs XUTs)
	GFFS=(RP_137_genes_TSS_Xu_2009.gff SAGA_TSS_Xu_2009_ORF_Ts_V64.gff TFIID_TSS_Xu_2009_ORF_Ts_V64.gff Xu_2009_CUTs_TSS_ONLY_V64.gff Xu_2009_SUTs_TSS_ONLY_V64.gff van_Dijk_2011_XUTs_V64_TSS_ONLY.gff)

	for j in `seq 0 $((${#GENE_CLASSES[@]}-1))`
	do
		GENE_CLASS=${GENE_CLASSES[$j]}
		GFF=${GFFS[$j]}

		OUT_DIR=$FACTOR_DIR/$GENE_CLASS

		if [ ! -d $OUT_DIR ]
			then
				python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR $NORM_DIR ../shared_files/$GFF
		fi
		python ../scripts/composite_plots.py -w 20 $OUT_DIR
	done

done
