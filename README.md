# vinayachandran_2017_figs

# Prerequisites

* linux
* python >= 2.7
    * numpy
    * scipy
    * matplotlib
    * scikit-learn
* java
* perl
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)
* [MEME suite](http://meme-suite.org/doc/install.html?man_type=web)

# Use
Each directory contains shell scripts to create each figure panel. For example, go to Fig4 and run fig4a.sh to create Fig. 4A. 

Most images are automatically created. However, heatmaps must be created manually. Open a sorted CDT (see below) with [TreeView](http://jtreeview.sourceforge.net/). Settings > Pixel Settings. Set global X and Y to Fill. Set zero to white. Adjust contrast as needed. 

Note: the shell scripts check whether output files exist before running certain steps. So be sure to remove bad output and empty directories before re-running a script.

# Output

script | output
--- | ---
fig1.sh | gene_CDT/factor/_composite/composite_plots_all_factors.svg; e.g. HSP42_CDT/Htz1/_composite/composite_plot_all_factors.svg
fig2a.sh | a_CDT/sorted/IDsacCer3_.cdt; e.g. a_CDT/sorted/50519sacCer3_.cdt
fig2b.sh | b_CDT/sorted/IDsacCer3_.cdt; e.g. b_CDT/sorted/50417sacCer3_.cdt
fig3.sh | CDT/sorted/IDsacCer3_.cdt; e.g. CDT/sorted/50428sacCer3_.cdt
fig4a.sh | a_divergent_CDT/sorted/IDsacCer3_.cdt; a_non-divergent_CDT/sorted/IDsacCer3_.cdt; e.g. a_divergent_CDT/sorted/59806sacCer3_.cdt
fig4b.sh | b_CDT/_composite/composite_plot_all_factors.svg
fig4c.sh | Fig4C.png
fig5a.sh | a_CDT/_composite/composite_plot_all_factors.svg
fig5b.sh | b_CDT/_composite/composite_plot_all_factors.svg
fig6a.sh | a_CDT/_composite/composite_plot_all_factors.svg
fig6b.sh | b_CDT/_composite/composite_plot_all_factors.svg
fig7a.sh | a_CDT/_composite/composite_plot_all_factors.svg
fig7b.sh | b_CDT/sorted/IDsacCer3_.cdt; e.g. b_CDT/sorted/51831sacCer3_.cdt
sup1a.sh | Sup1a.png
sup1b.sh | IDsacCer3_.png; e.g. 50428sacCer3_.png
sup2.sh | tab_files/factor/geneclass/_composite/composite_plot_all_factors.svg; e.g. tab_files/Spt3/RP/_composite/composite_plot_all_factors.svg
sup3a.sh | tab_files_a/factor/gene_classification/_composite/composite_plot_all_factors.svg; e.g. tab_files_a/Hsf1/activated/_composite/composite_plot_all_factors.svg
sup3b.sh | b_CDT/sorted/IDsacCer3_.cdt; e.g. b_CDT/sorted/50428sacCer3_.cdt
sup4.sh | bottom25/_composite/composite_plot_all_factors.svg; top25/_composite/composite_plot_all_factors.svg
sup8.sh | CDT/sorted/IDsacCer3_.cdt; e.g. CDT/sorted/53811sacCer3_.cdt
