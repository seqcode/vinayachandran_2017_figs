set -e

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

GFF=../shared_files/Xu_2009_ORF_TSS_TES_V64.gff

if [ ! -e $GFF ]
	then
		wget 	#TODO: add url
		mv Xu_2009_ORF_TSS_TES_V64.gff.gz ../shared_files
		gunzip ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff.gz
fi

mkdir -p stranded_tab_files
mkdir -p unstranded_tab_files

STRANDED_IDS=(50525 50526 53301 53302 50428 50429 50519 50520 53407 53408 53824 53825 53814 53815 53811 53812 53809 53810)
UNSTRANDED_IDS=(56422 56423)

for STRANDED_ID in "${STRANDED_IDS[@]}"
do
	cp $ALL_TAB/$STRANDED_ID"sacCer3".tab stranded_tab_files
done

for UNSTRANDED_ID in "${UNSTRANDED_IDS[@]}"
do
	cp $ALL_TAB/$UNSTRANDED_ID"sacCer3".tab unstranded_tab_files
done

if [ ! -d stranded_tab_files/Normalized_tab_files ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py --stranded stranded_tab_files ../shared_files/sacCer3.chrom.sizes
fi

if [ ! -d unstranded_tab_files/Normalized_tab_files ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py unstranded_tab_files ../shared_files/sacCer3.chrom.sizes
fi

GENES=(HSP42 RPL3 REB1)
XU_IDS=(YDR171W YOR063W YBR049C)

STRANDED_YLIMS=(20 60 150 120 20 30 40 40 80)
UNSTRANDED_YLIMS=(5 5 500)

for i in `seq 0 $((${#GENES[@]}-1))`
do
	GENE=${GENES[$i]}
	XU_ID=${XU_IDS[$i]}
	
	if [ ! -e ../shared_files/$GENE.gff ]
		then
			cat $GFF | awk -v var=$XU_ID '$9 == var {print $0}' > ../shared_files/$GENE.gff
	fi

	for j in `seq 0 $((${#STRANDED_IDS[@]}-1))`
	do
		ID=${STRANDED_IDS[$j]}
		YLIM=${STRANDED_YLIMS[$j]}
		mkdir -p $ID
		cp stranded_tab_files/Normalized_tab_files/$ID"_forward.tab" $ID
		cp stranded_tab_files/Normalized_tab_files/$ID"_reverse.tab" $ID


		if [ ! -d $ID/$GENE"_CDT" ]
			then
				mkdir $ID/$GENE"_CDT"
				python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $ID/$GENE"_CDT" $ID ../shared_files/$GENE.gff
		fi

		python ../scripts/Composite_for_many_factors_one_plot.py -w 20 -y $YLIM $ID/$GENE"_CDT"

	done

	for j in `seq 0 $((${#UNSTRANDED_IDS[@]}-1))`
	do
		ID=${UNSTRANDED_IDS[$j]}
		YLIM=${UNSTRANDED_YLIMS[$j]}
		mkdir -p $ID
		mv unstranded_tab_files/Normalized_tab_files/$ID"_forward.tab" $ID
		mv unstranded_tab_files/Normalized_tab_files/$ID"_reverse.tab" $ID


		if [ ! -d $ID/$GENE"_CDT" ]
			then
				mkdir $ID/$GENE"_CDT"
				python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $ID/$GENE"_CDT" $ID ../shared_files/$GENE.gff
		fi

		python ../scripts/composite_for_many_factors_one_plot.py -w 20 -y $YLIM $ID/$GENE"_CDT"

	done
done 
