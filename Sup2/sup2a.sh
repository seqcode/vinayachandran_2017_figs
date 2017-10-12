set -e

NORM_DIR=tab_files_a/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_a ../shared_files/sacCer3.chrom.sizes
fi

#OUT_DIRS=(RP SAGA TFIID CUTs SUTs XUTs)
#GFFS=(RP_137_genes_TSS_Xu_2009.gff SAGA_TSS_Xu_2009_ORF_Ts_V64.gff TFIID_TSS_Xu_2009_ORF_Ts_V64.gff Xu_2009_CUTs_TSS_ONLY_V64.gff Xu_2009_SUTs_TSS_ONLY_V64.gff van_Dijk_2011_XUTs_V64_TSS_ONLY.gff)
OUT_DIRS=(SAGA RP)		#TEST
GFFS=(SAGA_TSS_Xu_2009_ORF_Ts_V64.gff RP_137_genes_TSS_Xu_2009.gff)

for i in `seq 0 $((${#OUT_DIRS[@]}-1))`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR tab_files_a/Normalized_tab_files ../shared_files/$GFF
	fi
	python ../scripts/composite_plots.py -w 20 --shaded $OUT_DIR
done
