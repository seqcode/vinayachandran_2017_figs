set -e

NORM_DIR=tab_files_c/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_c ../shared_files/sacCer3.chrom.sizes
fi

#OUT_DIRS=(RP SAGA_induced SAGA_repressed SAGA_no_change TFIID_induced TFIID_repressed TFIID_no_change)
#GFFS=(RP_137_genes_TSS_Xu_2009_PURE_SET.gff SAGA-activated_TSS_Xu_2009_ORF_PURE_SET.gff SAGA-repressed_TSS_Xu_2009_ORF_PURE_SET.gff SAGA-nochange_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-activated_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-nochange_TSS_Xu_2009_ORF_PURE_SET.gff TFIID-repressed_TSS_Xu_2009_ORF_PURE_SET.gff)
OUT_DIRS=(SAGA_induced)	#TEST
GFFS=(SAGA-activated_TSS_Xu_2009_ORF_PURE_SET.gff)

for i in `seq 0 $((${#OUT_DIRS[@]}-1))`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 -o $OUT_DIR $NORM_DIR ../shared_files/$GFF
	fi

	python ../scripts/calculate_sum_from_columns.py .tab 50 50
done
