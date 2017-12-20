set -e

sh ../scripts/get_chrom_sizes.sh

TAB_DIR=tab_files_c
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

CDT_DIR=c_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 100 -d 100 -o $CDT_DIR $TAB_DIR ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff
fi

TRANSCRIPTION_FILE=yassour.tsv

if [ ! -e $TRANSCRIPTION_FILE ] 
	then
		wget http://compbio.cs.huji.ac.il/RNASeq/SuppSite/Supp_Tables_files/SuppTable3_Genes.xls.gz
		gunzip SuppTable3_Genes.xls.gz
		cat SuppTable3_Genes.xls | awk -F "\t" '$1 != "GENE" {print $1"\t"$8}' > $TRANSCRIPTION_FILE
		rm SuppTable3_Genes.xls
fi

python correlate_with_occupancy.py $CDT_DIR $TRANSCRIPTION_FILE "Transcription (mRNA)"
