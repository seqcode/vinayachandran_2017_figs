set -e

WD=$PWD

SHARED_FILES=../shared_files
if [ ! -d $SHARED_FILES ]
	then
		mkdir $SHARED_FILES
fi

CHROM_INFO=../shared_files/sacCer3.chrom.sizes
if [ ! -e $CHROM_INFO ]
	then
		cd $SHARED_FILES
		wget http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
		cd $WD
fi

python ../scripts/quantile_norm_singlebase_bin.py tab_files_a $CHROM_INFO

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

	python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR tab_files_a/Normalized_tab_files $GFF
	python ../scripts/composite_plots.py -w 20 --shaded $OUT_DIR
done
