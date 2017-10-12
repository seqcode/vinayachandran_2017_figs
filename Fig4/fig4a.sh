set -e

NORM_DIR=tab_files_a/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_a ../shared_files/sacCer3.chrom.sizes
fi

OUT_DIRS=(a_divergent_CDT a_non-divergent_CDT)
GFFS=(Hsf1-union-Xu-TSS-divergent-upstream.gff Hsf1-union-Xu-TSS-sortby-distance.gff)

for i in `seq 0 $((${#OUT_DIRS[@]}-1))`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $OUT_DIR $NORM_DIR ../shared_files/$GFF
	python ../scripts/sort_cdt_by_given_file.py -o 2 $OUT_DIR ../shared_files/$GFF
done
