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
		mkdir control	#negative control is in separate directory so that it won't be processed
		mv GSM2865467_15941sacCer3.tab control
		cd ..
fi

#call peaks
python ../scripts/chipexo/genetrack/genetrack.py -s 5 -e 20 tab_files_b/GSM2601075_53302sacCer3.tab > genetrack_peaks.gff

#filter out singletons
python filter_singletons.py genetrack_peaks.gff filtered.gff

#pair peaks
python cwpair2.py -u 80 -d 80 -b 2 filtered.gff

#get top 500
sort --reverse -k 6 -n cwpair_output_mode_f0u80d80b2/S_filtered.gff | head -500 > top_500.gff

#get 80bp windows
bedtools slop -i top_500.gff -g ../shared_files/sacCer3.chrom.sizes -b 40 > top_500_80bp.gff

#get sacCer3 fasta
if [ ! -e twoBitToFa ]
	then
		wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
		chmod +x twoBitToFa
fi

if [ ! -e sacCer3.fa ]
	then
		./twoBitToFa http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit stdout > sacCer3.fa
		python ../scripts/chipexo/genetrack/chrtrans.py sacCer3.fa
		mv roman_to_numeric/sacCer3.fa sacCer3.fa
		rm -r roman_to_numeric
fi

#convert to fasta
bedtools getfasta -fi sacCer3.fa -bed top_500_80bp.gff > top_500_80bp.fa

#run MEME
meme top_500_80bp.fa -dna -revcomp

#find all genomic sites with motif
fimo --thresh 0.0001 meme_out/meme.txt sacCer3.fa

#center binding locations on significant motifs
perl matchGFFwithFIMO.pl cwpair_output_mode_f0u80d80b2/S_filtered.gff fimo_out/fimo.gff

#calculate enrichment over control
java -jar SignificanceTester_pugh_java1.7.jar --geninfo ../shared_files/sacCer3.chrom.sizes --format IDX --expt tab_files_b/GSM2601075_53302sacCer3.tab --ctrl tab_files_b/control/GSM2865467_15941sacCer3.tab --gff cwpair_output_mode_f0u80d80b2/S_filtered_withmotif.gff --q 0.05 --minfold 1

NORM_DIR=$TAB_DIR/Normalized_tab_files

if [ ! -d $NORM_DIR ]
	then
		python ../scripts/quantile_norm_singlebase_bin.py $TAB_DIR ../shared_files/sacCer3.chrom.sizes
fi

CDT_DIR=b_CDT

if [ ! -d $CDT_DIR ]
	then
		python ../scripts/map_shifted_tags_to_ref.py -u 400 -d 400 -o $CDT_DIR $NORM_DIR cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold1.0/S_filtered_withmotif_signif_EXPERIMENT.gff
fi

python ../scripts/composite_plots.py -w 20 --normalize $CDT_DIR
