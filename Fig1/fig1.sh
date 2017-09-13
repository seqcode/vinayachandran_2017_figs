set -e

python ../scripts/quantile_norm_singlebase_bin_strandness.py . 
python ../scripts/Composite_for_many_factors_one_plot.py -w 20 . 
