set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

OUT_DIRS=(activated nochange repressed_rp)
GFFS=(activated.gff nochange.gff repressed_rp.gff)

for i in `seq 0 ${#OUT_DIRS[@]}`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $OUT_DIR tab_files/Normalized_tab_files $GFF
	python ../scripts/Composite_plots_Vinesh\(1\).py -w 20 $OUT_DIR
done
