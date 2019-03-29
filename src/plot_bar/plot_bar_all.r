
library(data.table)

plot_figures <- function(case, output_name) {

fc_threshold <- sort(abs(case$log2FC), decreasing=T)[50]
case_filtered <- case[which(abs(case$log2FC) >= fc_threshold), ]

case_sorted <- case_filtered[order(case_filtered$log2FC), ]

case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene <- paste0( "***   ", case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene)

my_col <- rep("gray", dim(case_sorted)[2])
my_col[which(case_sorted$log2FC > 1)] <- "orange"
my_col[which(case_sorted$log2FC < -1)] <- "blue"

my_xmax <- round(max(abs(case_sorted$log2fc), na.rm=T)) + 3


pdf(paste0("barplot_", output_name, "_all_top50.pdf"))

barplot(case_sorted$log2FC, names.arg=case_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2FC", xlim=c(-my_xmax, my_xmax), cex.names=0.5, main=paste0(dim(case_filtered)[1], " proteins"))

abline(v=fc_threshold, lty=3)
abline(v=-fc_threshold, lty=3)
dev.off()

}


list_tsv <- system(command="ls *_preprocessed.tsv", intern=T)
cosmic <- fread("~/project/phrt/TpViz/docker-tp-viz/data/Cancer_Census_all_072018_COSMIC.csv")

for(i in 1:length(list_tsv)) {

filename <- gsub("_preprocessed.tsv", "", list_tsv[i])

raw <- fread(list_tsv[i], dec=",")

plot_figures(raw, filename)


#plot_figures(raw)

}
