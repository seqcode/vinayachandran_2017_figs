set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

OUT_DIRS=(RP SAGA TFIID CUTs SUTs XUTs)
GFFS=(RP_137_genes_TSS_Xu_2009.gff SAGA_TSS_Xu_2009_ORF_Ts_V64.gff TFIID_TSS_Xu_2009_ORF_Ts_V64.gff Xu_2009_CUTs_TSS_ONLY_V64.gff Xu_2009_SUTs_TSS_ONLY_V64.gff van_Dijk_2011_XUTs_V64_TSS_ONLY.gff)

for i in `seq 0 ${#OUT_DIRS[@]}`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR tab_files/Normalized_tab_files $GFF
	python ../scripts/composite_plots_shaded.py -w 20 $OUT_DIR
done
