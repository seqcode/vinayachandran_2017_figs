set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files_a 

OUT_DIRS=(divergent non-divergent)
GFFS=(Hsf1-union-Xu-TSS-divergent-upstream.gff Hsf1-union-Xu-TSS-sortby-distance.gff)

for i in `seq 0 ${#OUT_DIRS[@]}`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $OUT_DIR tab_files_a/Normalized_tab_files $GFF
	python ../scripts/sort_cdt_by_given_file.py -o 2 $OUT_DIR Hsf1-union-Xu-TSS-sortby-distance.gff
done
