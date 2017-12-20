set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files_b
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

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 100 -d 100 -o $CDT_DIR $NORM_DIR ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff
fi

TRANSCRIPTION_RATE_FILE=holstege.tsv

if [ ! -e $TRANSCRIPTION_RATE_FILE ] 
	then
		wget http://younglab.wi.mit.edu/pub/data/orf_transcriptome.txt
		cat orf_transcriptome.txt | awk '$1 != "ORF" && $4 != "#N/A" {print $1"\t"$4}' > $TRANSCRIPTION_RATE_FILE
		rm orf_transcriptome.txt
fi

python correlate_with_occupancy.py $CDT_DIR $TRANSCRIPTION_RATE_FILE "Transcription rate (mRNA per hr)"
