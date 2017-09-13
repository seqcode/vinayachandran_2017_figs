set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 tab_files/Normalized_tab_files $GFF	#TODO: where is gff
fi

python ../scripts/extract_tag_occupancy.py $FN 100 300	#TODO: what is fn
