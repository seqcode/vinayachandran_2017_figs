set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 

CDT_DIR=_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 200 -d 200 tab_files/Normalized_tab_files SAGA_TSS_Xu_2009_ORF_Ts_V64.gff

fi

python ../scripts/composite_plots_shaded.py -w 21 $CDT_DIR