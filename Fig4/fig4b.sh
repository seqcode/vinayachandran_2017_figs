set -e

NORM_DIR=tab_files_b/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		#python ../scripts/quantile_norm_singlebase_bin.py tab_files_b ../shared_files/sacCer3.chrom.sizes
		python ~/Downloads/quantile_norm_singlebase_bin.py -g ../shared_files/sacCer3.chrom.sizes tab_files_b 	#TEST
fi

#CDT_DIR=b_CDT
CDT_DIR=tab_files_b/_CDT	#TEST
GFF=../shared_files/Hsf1-union-midpoint

if [ ! -e $GFF.gff ]
	then
		sh ../scripts/csv_to_gff.sh $GFF
fi

if [ ! -d $CDT_DIR ]
	then
		#python ../scripts/map_shifted_tags_to_ref.py -u 400 -d 400 -o $CDT_DIR tab_files_b/Normalized_tab_files $GFF.gff
		python ~/Downloads/map_shifted_tags_to_ref.py -u 400 -d 400 -r $GFF.gff -l /usr/bin/ tab_files_b/Normalized_tab_files	#TEST
fi

#python ../scripts/composite_plots.py -w 20 --shaded --normalize $CDT_DIR
python ~/Downloads/Composite_plots_shaded_bxp_Vinesh.py -w 20 $CDT_DIR	#TEST
