set -e

if [ ! -d ../shared_files ]
	then
		mkdir ../shared_files
fi

if [ ! -d ../scripts/chipexo ]
	then
		git clone git@github.com:ialbert/chipexo.git
		mv chipexo ../scripts
fi

if [ ! -e ../shared_files/sacCer3.chrom.sizes ]
	then
		wget http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
		python ../scripts/chipexo/genetrack/chrtrans.py sacCer3.chrom.sizes
		mv roman_to_numeric/sacCer3.chrom.sizes ../shared_files/sacCer3.chrom.sizes
		rm -r roman_to_numeric
fi
