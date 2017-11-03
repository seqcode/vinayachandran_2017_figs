set -e

if [ ! -d tab_files_c/Normalized_tab_files ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files_c ../shared_files/sacCer3.chrom.sizes
fi

python fig4c.py
