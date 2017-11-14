set -e

sh ../scripts/get_chrom_sizes.sh

ALL_TAB=../GSE98573_RAW

if [ ! -d $ALL_TAB ]
	then 
		tar xvf $ALL_TAB.tar
fi

IDS=(53301 59801 59804 50428 53302 59802 59805 50429 53303 59803 59806 50430 53824 53825 53826)
TAB_DIR=tab_files_a

if [ ! -d $TAB_DIR ]
	then
		mkdir $TAB_DIR

		for ID in "${IDS[@]}"
		do
			cp $ALL_TAB/$ID"sacCer3".tab $TAB_DIR
		done
fi

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

OUT_DIRS=(a_divergent_CDT a_non-divergent_CDT)
GFFS=(Hsf1-union-Xu-TSS-divergent-upstream.gff Hsf1-union-Xu-TSS-sortby-distance.gff)

for i in `seq 0 $((${#OUT_DIRS[@]}-1))`
do
	OUT_DIR=${OUT_DIRS[$i]}
	GFF=${GFFS[$i]}
	if [ ! -d $OUT_DIR ]
		then
			mkdir $OUT_DIR
	fi

	if [ ! -e $GFF ]
		then
			wget 	#TODO: add url
			gunzip $GFF.gz
	fi

	python ../scripts/map_shifted_tags_to_ref.py -u 1000 -d 1000 -o $OUT_DIR $NORM_DIR $GFF
	python ../scripts/sort_cdt_by_given_file.py -o 2 $OUT_DIR $GFF
done
