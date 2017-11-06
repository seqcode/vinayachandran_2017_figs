set -e

FACTORS=(Ssl2 Rpb3)
YLIMS=(2000 400)
CONDITIONS=(activated nochange repressed_rp)

for i in `seq 0 $((${#FACTORS[@]}-1))`
do
	FACTOR=${FACTORS[$i]}
	YLIM=${YLIMS[$i]}
	TAB_DIR=tab_files_a/$FACTOR
	NORM_DIR=$TAB_DIR/Normalized_tab_files
	if [ ! -d $NORM_DIR ]
		then
			python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
	fi

	for CONDITION in "${CONDITIONS[@]}"
	do
		CDT_DIR=$NORM_DIR/$CONDITION
		if [ ! -d $CDT_DIR ]
			then
				python ../scripts/map_shifted_tags_to_ref.py -u 500 -d 500 -o $CDT_DIR $NORM_DIR ../shared_files/$CONDITION.gff
		fi

		#python ../scripts/composite_plots.py -w 20 -y $YLIM $CDT_DIR
		python ../scripts/composite_plots.py -w 20 $CDT_DIR
	done
done
