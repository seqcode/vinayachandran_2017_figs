set -e

python ../scripts/quantile_norm_singlebase_bin.py tab_files 
python ../scripts/map_shifted_tags_to_ref.py -u 2000 -d 2000 tab_files/Normalized_tab_files 
python ../scripts/sort_CDT_by_given_file.py -o 2 tab_files/Normalized_tab_files/_CDT
