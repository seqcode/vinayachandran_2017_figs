set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files_a
TAR=GSE98573_RAW.tar

if [ ! -d $TAB_DIR ]
	then 
		mkdir $TAB_DIR	
		mv $TAR $TAB_DIR
		cd $TAB_DIR
		tar xvf $TAR
		rm $TAR
		gunzip *.gz
		cd ..
fi

FACTORS=(Spt3 Spt15 Sua7 Taf1 Ssl2 Rpb3 Hsf1)
YLIMS=(3000 2500 2000 900 2000 400 1600)
ALL_IDS=(50525,50526,50527,50528,50529,50530 53319,53320,53321 50428,50429,50430 50531,50532,50533,50534,50535,50536 50519,50520,50521,50522,50523,50524 51835,51836,51837,51838,51839,51840 59804,59805,59806)
CONDITIONS=(activated nochange repressed_rp)

for i in `seq 0 $((${#FACTORS[@]}-1))`
do
	FACTOR=${FACTORS[$i]}
	YLIM=${YLIMS[$i]}
	FACTOR_DIR=$TAB_DIR/$FACTOR
	NORM_DIR=$FACTOR_DIR/Normalized_tab_files
	IDS=${ALL_IDS[$i]}
	IDS_ARRAY=(${IDS//,/ })
	
	if [ ! -d $FACTOR_DIR ]
		then
			mkdir -p $FACTOR_DIR
			for ID in "${IDS_ARRAY[@]}"
			do
				mv $TAB_DIR/*_$ID"sacCer3".tab $FACTOR_DIR
			done
	fi

	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py $FACTOR_DIR ../shared_files/sacCer3.chrom.sizes
	fi

	for j in `seq 0 $((${#CONDITIONS[@]}-1))`
	do
		CONDITION=${CONDITIONS[$j]}

		CDT_DIR=$FACTOR_DIR/$CONDITION
		if [ ! -d $CDT_DIR ]
			then
				python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR ../shared_files/$CONDITION.gff
		fi

		python ../scripts/composite_plots.py -w 20 -y $YLIM $CDT_DIR
	done
done
