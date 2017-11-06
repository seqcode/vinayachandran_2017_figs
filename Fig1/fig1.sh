set -e

if [ ! -d stranded_tab_files/Normalized_tab_files ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py --stranded stranded_tab_files ../shared_files/sacCer3.chrom.sizes
		#python ~/Downloads/quantile_norm_singlebase_bin_strandness.py -g ../shared_files/sacCer3.chrom.sizes stranded_tab_files
fi

#if [ ! -d unstranded_tab_files/Normalized_tab_files ]
#	then
#		python ../scripts/quantile_norm_singlebase_bin.py unstranded_tab_files ../shared_files/sacCer3.chrom.sizes
#fi

#GENES=(HSP42 RPL3 REB1)
GENES=(REB1)
#XU_IDS=(YDR171W YOR063W YBR049C)
XU_IDS=(YBR049C)

YLIM=40	#TODO: don't hardcode

#for ID in 53811sacCer3 53812sacCer3
for ID in 53811sacCer3
do
	mkdir $ID
	#mv stranded_tab_files/Normalized_tab_files/$ID"_forward.tab" $ID
	#mv stranded_tab_files/Normalized_tab_files/$ID"_reverse.tab" $ID
	cp stranded_tab_files/Normalized_tab_files/$ID"_forward.tab" $ID
	cp stranded_tab_files/Normalized_tab_files/$ID"_reverse.tab" $ID
	for i in `seq 0 $((${#GENES[@]}-1))`
	do
		GENE=${GENES[$i]}
		XU_ID=${XU_IDS[$i]}
		
		if [ ! -e ../shared_files/$GENE.gff ]
			then
				cat ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff | awk -v var=$XU_ID '$9 == var {print $0}' > ../shared_files/$GENE.gff
		fi

		if [ ! -d $ID/$GENE"_CDT" ]
			then
				mkdir $ID/$GENE"_CDT"
				python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $ID/$GENE"_CDT" $ID ../shared_files/$GENE.gff
				#python ~/Documents/vinesh/map_to_Reb1/map_shifted_tags_to_ref.py -u 1000 -d 1000 -r ../shared_files/$GENE.gff -l /usr/bin/ $ID 
		fi

		python ../scripts/Composite_for_many_factors_one_plot.py -w 20 -y $YLIM $ID/$GENE"_CDT"

	done
done

#rm -r stranded_tab_files/Normalized_tab_files
