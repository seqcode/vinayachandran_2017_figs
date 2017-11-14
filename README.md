# vinayachandran_2017_figs

# Prerequisites

* python >= 2.7
    * numpy
    * scipy
    * matplotlib
    * scikit-learn
* java
* perl
* bedtools (http://bedtools.readthedocs.io/en/latest/content/installation.html)
* MEME suite (http://meme-suite.org/doc/install.html?man_type=web)

# Use
Each directory contains shell scripts to create each figure panel. For example, run Fig4/fig4a.sh to create Fig. 4A. 

Most images are automatically created. However, heatmaps must be created manually. Open a sorted CDT (see below) with TreeView (http://jtreeview.sourceforge.net/). Settings > Pixel Settings. Set global X and Y to Fill. Set zero to white. Adjust contrast as needed. 

Note: the shell scripts check whether output files exist before running certain steps. So be sure to remove bad output and empty directories before re-running the script.

# Output

script | output
--- | ---
fig1.sh | gene_CDT/factor/_composite/composite_plots_all_factors.svg; e.g. HSP42_CDT/Htz1/_composite/composite_plot_all_factors.svg
