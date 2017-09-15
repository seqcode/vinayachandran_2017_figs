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

OUT_DIRS=(RP SAGA_induced SAGA_repressed SAGA_no_change TFIID_induced TFIID_repressed TFIID_no_change)
GFFS=(RP_137_genes_TSS_Xu_2009_PURE_SET.gff SAGA-activated_TSS_Xu_2009_ORF_PURE_SET.gff SAGA-nochange_TSS_Xu_2009_ORF_PURE_SET.gff SAGA-repressed_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-activated_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-nochange_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-repressed_TSS_Xu_2009_ORF_PURE_SET.gff)

for i in `seq 0 ${#OUT_DIRS[@]}`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 -o $OUT_DIR tab_files_c/Normalized_tab_files $GFF
	python ../scripts/calculate_sum_from_columns.py -s .tab -d 50 -u 50
done
