set -e

for tab_file in tab_files_a/*
do
	Rscript Ssl2_DEseq_replicates.R $tab_file
	python calculate_sum_from_columns.py .tab
done
