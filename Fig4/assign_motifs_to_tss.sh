set -e

#create motif sequences file
bedtools getfasta -fi sacCer3.fa -bed cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff > motifs.fa

#sort files
bedtools sort -i cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff > peaks_with_motif_sorted.gff
bedtools sort -i ../shared_files/Xu_2009_ORF_TSS_TES_V64.gff > ../shared_files/Xu_2009_ORF_TSS_TES_V64_sorted.gff

#find closest TSS within 500 bp for each motif
bedtools closest -a peaks_with_motif_sorted.gff -b ../shared_files/Xu_2009_ORF_TSS_TES_V64_sorted.gff -d > closest_tss.gff
cat closest_tss.gff | awk '$19 < 500 {print $0}' > closest_tss_within_500bp.gff

#create assignment list
python assign_motifs_to_tss.py

#create GFF of non-duplicate TSSs
cat assignment_list.tsv | awk '{print $1}' > tss_nodups.txt
cat closest_tss_within_500bp.gff | awk '{print $10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}' > tss_only.gff
python ../Sup4/genes_from_gff.py tss_only.gff tss_nodups.txt
bedtools sort -i tss_nodups.gff > tss_nodups_sorted.gff 

#create divergent list
bedtools closest -a tss_nodups_sorted.gff -b ../shared_files/Xu_2009_ORF_TSS_TES_V64_sorted.gff -S > opposite_strand.gff
python check_for_divergent.py

#create list of all TSSs
cat divergent_list.tsv | awk '$2 != "None" {print $2}' > divergent_names.txt
python ../Sup4/genes_from_gff.py ../shared_files/Xu_2009_ORF_TSS_TES_V64_sorted.gff divergent_names.txt 
cp tss_nodups.gff all_tss.gff
cat divergent.gff >> all_tss.gff

cat tss_nodups.gff | awk '{print $9}' > tss_nodups_names.txt
cp tss_nodups_names.txt all_tss_names.txt
cat divergent_names.txt >> all_tss_names.txt

TFIIB_FLANK=150
HSF1_FLANK=100
SPT3_FLANK=150

#TFIIB MHS occupancy
python ../scripts/extract_tag_occupancy.py tab_files_a/50428sacCer3.tab all_tss.gff ../shared_files/sacCer3.chrom.sizes TFIIB_MHS_occupancy.txt $TFIIB_FLANK $TFIIB_FLANK
paste all_tss_names.txt TFIIB_MHS_occupancy.txt > TFIIB_MHS_occupancy.tsv

#TFIIB HS3 occupancy
python ../scripts/extract_tag_occupancy.py tab_files_a/50429sacCer3.tab all_tss.gff ../shared_files/sacCer3.chrom.sizes TFIIB_HS3_occupancy.txt $TFIIB_FLANK $TFIIB_FLANK
paste all_tss_names.txt TFIIB_HS3_occupancy.txt > TFIIB_HS3_occupancy.tsv

#create motif IDs
cat cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff | awk '{print $1":"$4"-"$5}' > motif_ids.txt

#Hsf1 MHS occupancy
python ../scripts/extract_tag_occupancy.py tab_files_a/53301sacCer3.tab cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff ../shared_files/sacCer3.chrom.sizes Hsf1_MHS_occupancy.txt $HSF1_FLANK $HSF1_FLANK
paste motif_ids.txt Hsf1_MHS_occupancy.txt > Hsf1_MHS_occupancy.tsv
 
#Hsf1 HS3 occupancy
python ../scripts/extract_tag_occupancy.py tab_files_a/53302sacCer3.tab cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff ../shared_files/sacCer3.chrom.sizes Hsf1_HS3_occupancy.txt $HSF1_FLANK $HSF1_FLANK
paste motif_ids.txt Hsf1_HS3_occupancy.txt > Hsf1_HS3_occupancy.tsv

#Spt3 MHS occupancy
python ../scripts/extract_tag_occupancy.py ../Sup2/tab_files/Spt3/50525sacCer3.tab cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff ../shared_files/sacCer3.chrom.sizes Spt3_MHS_occupancy.txt $SPT3_FLANK $SPT3_FLANK
paste motif_ids.txt Spt3_MHS_occupancy.txt > Spt3_MHS_occupancy.tsv

#Spt3 HS3 occupancy
python ../scripts/extract_tag_occupancy.py tab_files_b/50526sacCer3.tab cwpair_output_mode_f0u80d80b2/signif_w50_q5.00e-02_minfold2.0/S_filtered_withmotif_signif_EXPERIMENT.gff ../shared_files/sacCer3.chrom.sizes Spt3_HS3_occupancy.txt $SPT3_FLANK $SPT3_FLANK
paste motif_ids.txt Spt3_HS3_occupancy.txt > Spt3_HS3_occupancy.tsv

python fill_table_s1.py
