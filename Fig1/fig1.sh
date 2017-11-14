set -e

sh ../scripts/get_chrom_sizes.sh

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

GFF=Xu_2009_ORF_TSS_TES_V64.gff

if [ ! -e ../shared_files/$GFF ]
	then
		wget 	#TODO: add url
		mv $GFF.gz ../shared_files
		gunzip ../shared_files/$GFF.gz
fi

TAB_DIR=tab_files

STRANDED_FACTORS=(SAGA Hsf1 TFIIB TFIIH FACT PolII Ser7p Ser5p Ser2p)
ALL_STRANDED_IDS=([50525 50526] [53301 53302] [50428 50429] [50519 50520] [53407 53408] [53824 53825] [53814 53815] [53811 53812] [53809 53810])

for i in `seq 0 $((${#STRANDED_FACTORS[@]}-1))`
do
	FACTOR=${STRANDED_FACTORS[$i]}
	IDS=${ALL_STRANDED_IDS[$i]}

	FACTOR_DIR=$TAB_DIR/$FACTOR

	if [ ! -d $FACTOR_DIR ]
		then
			mkdir $FACTOR_DIR

			for ID in "${IDS[@]}"
			do
				cp $ALL_TAB/$ID"sacCer3".tab $FACTOR_DIR
			done
	fi

	NORM_DIR=$FACTOR_DIR/Normalized_tab_files
	
	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py --stranded $FACTOR_DIR ../shared_files/sacCer3.chrom.sizes
	fi

done

UNSTRANDED_FACTORS=(Htz1 PIP-seq)
ALL_UNSTRANDED_IDS=([50416 50417] [56422 56423])

for i in `seq 0 $((${#UNSTRANDED_FACTORS[@]}-1))`
do
	FACTOR=${UNSTRANDED_FACTORS[$i]}
	IDS=${ALL_UNSTRANDED_IDS[$i]}

	FACTOR_DIR=$TAB_DIR/$FACTOR

	if [ ! -d $FACTOR_DIR ]
		then
			mkdir $FACTOR_DIR

			for ID in "${IDS[@]}"
			do
				cp $ALL_TAB/$ID"sacCer3".tab $FACTOR_DIR
			done
	fi

	NORM_DIR=$FACTOR_DIR/Normalized_tab_files
	
	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py $FACTOR_DIR ../shared_files/sacCer3.chrom.sizes
	fi
done

GENES=(HSP42 RPL3 REB1)
XU_IDS=(YDR171W YOR063W YBR049C)

STRANDED_YLIMS=(20 60 150 120 20 30 40 40 80)
UNSTRANDED_YLIMS=(5 500)

for i in `seq 0 $((${#GENES[@]}-1))`
do
	GENE=${GENES[$i]}
	XU_ID=${XU_IDS[$i]}
	
	if [ ! -e ../shared_files/$GENE.gff ]
		then
			cat ../shared_files/$GFF | awk -v var=$XU_ID '$9 == var {print $0}' > ../shared_files/$GENE.gff
	fi

	for j in `seq 0 $((${#STRANDED_FACTORS[@]}-1))`
	do
		FACTOR=${STRANDED_FACTORS[$j]}
		YLIM=${STRANDED_YLIMS[$j]}
		FACTOR_DIR=$TAB_DIR/$FACTOR

		if [ ! -d $FACTOR_DIR/$GENE"_CDT" ]
			then
				mkdir $FACTOR_DIR/$GENE"_CDT"
				python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $FACTOR_DIR/$GENE"_CDT" $FACTOR_DIR/Normalized_tab_files ../shared_files/$GENE.gff
		fi

		python ../scripts/composite_plots.py -w 20 -y $YLIM --shaded $FACTOR_DIR/$GENE"_CDT"

	done

	for j in `seq 0 $((${#UNSTRANDED_FACTORS[@]}-1))`
	do
		FACTOR=${UNSTRANDED_FACTORS[$j]}
		YLIM=${UNSTRANDED_YLIMS[$j]}
		FACTOR_DIR=$TAB_DIR/$FACTOR

		if [ ! -d $FACTOR_DIR/$GENE"_CDT" ]
			then
				mkdir $FACTOR_DIR/$GENE"_CDT"
				python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $FACTOR_DIR/$GENE"_CDT" $FACTOR_DIR/Normalized_tab_files ../shared_files/$GENE.gff
		fi

		python ../scripts/composite_plots.py -w 20 -y $YLIM --shaded $FACTOR_DIR/$GENE"_CDT"

	done
done 
