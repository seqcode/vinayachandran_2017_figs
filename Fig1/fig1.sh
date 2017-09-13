set -e

python ../scripts/quantile_norm_singlebase_bin_strandness.py tab_files
python ../scripts/Composite_for_many_factors_one_plot.py -w 20 tab_files/Normalized_tab_files 
