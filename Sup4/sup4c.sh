set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

OUT_DIR=_CDT
if [ ! -d $OUT_DIR ]
	then
		mkdir $OUT_DIR
fi

python ../scripts/map_shifted_tags_to_ref.py -u 450 -d 450 -o $OUT_DIR tab_files/Normalized_tab_files $GFF	#TODO: where is gff
python ../scripts/composite_plots_shaded.py -w 20 $OUT_DIR
