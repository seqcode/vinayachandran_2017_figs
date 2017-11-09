set -e

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

TAB_DIR=tab_files_a

FACTORS=(Spt3 Spt15 Sua7 Taf1 Ssl2 Rpb3 Hsf1)
YLIMS=(3000 2500 2000 900 2000 400 1600)
ALL_IDS=([50525 50526 50527 50528 50529 50530] [53319 53320 53321] [50428 50429 50430] [50531 50532 50533 50534 50535 50536] [50519 50520 50521 50522 50523 50524] [51835 51836 51837 51838 51839 51840] [53301 59801 59804 53302 59802 59805 53303 59803 59806])
CONDITIONS=(activated nochange repressed_rp)
URLS=()	#TODO: add

for i in `seq 0 $((${#FACTORS[@]}-1))`
do
	FACTOR=${FACTORS[$i]}
	YLIM=${YLIMS[$i]}
	FACTOR_DIR=$TAB_DIR/$FACTOR
	NORM_DIR=$FACTOR_DIR/Normalized_tab_files
	IDS=${ALL_IDS[$i]}
	
	if [ ! -d $FACTOR_DIR ]
		then
			mkdir -p $FACTOR_DIR
			for ID in "${IDS[@]}"
			do
				cp $ALL_TAB/$ID"sacCer3".tab $FACTOR_DIR
			done
	fi

	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py $FACTOR_DIR ../shared_files/sacCer3.chrom.sizes
	fi

	for j in `seq 0 $((${#CONDITIONS[@]}-1))`
	do
		URL=${URLS[$j]}
		CONDITION=${CONDITIONS[$j]}

		GFF=$CONDITION.gff
		if [ ! -e ../shared_files/$GFF ]
			then
				wget $URL
				mv $GFF.gz ../shared_files
				gunzip ../shared_files/$GFF.gz
		fi

		CDT_DIR=$NORM_DIR/$CONDITION
		if [ ! -d $CDT_DIR ]
			then
				python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR ../shared_files/$GFF
		fi

		python ../scripts/composite_plots.py -w 20 -y $YLIM $CDT_DIR
	done
done
