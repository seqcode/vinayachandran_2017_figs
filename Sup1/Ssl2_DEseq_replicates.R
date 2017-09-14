## abc.txt is a matrix of row containing gene names and col containing sample names, each element in the matrix represents the sum of tags/reads in a given interval around TSS (+/- 100bp) for a given gene and sample. For PCA plot in the figure, TFIIH dataset was used for all time points of heat shock (0,3,6,9,12,15).

library(ggplot2)
dat = read.table('abc.txt',header = TRUE, sep="\t")
# Assuming that gene names are in column 1
dat = data.frame(dat, row.names=1)
pca = prcomp(t(dat), scale=TRUE, center=TRUE)
pcs = data.frame(PC1=pca$x[,"PC1"], PC2=pca$x[,"PC2"], treatment=colnames(dat))
ggplot(pcs, aes(x=PC1, y=PC2, color = treatment)) + geom_point(size = 3)
