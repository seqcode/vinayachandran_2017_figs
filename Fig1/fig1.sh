set -e

if [ ! -d tab_files ]
	then
		mkdir tab_files
		cat ../shared_files/Figure_Script_Table.csv | awk -F "," '$4 == 1 {print $11}' > tab_files_list
		python ../scripts/fetch_from_gpfs.py Fig1
fi

for f in tab_files/*.gz
do
	gunzip $f
done

if [ ! -d tab_files/Normalized_tab_files ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py tab_files ../shared_files/sacCer3.chrom.sizes --stranded
fi

if [ ! -d tab_files/_CDT ]
	then
		python ../scripts/map_shifted_tags_to_ref.py tab_files/Normalized_tab_files ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff
fi

python ../scripts/Composite_for_many_factors_one_plot.py -w 20 tab_files/Normalized_tab_files 
