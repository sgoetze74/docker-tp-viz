
library(data.table)

plot_figures <- function(case, output_name) {

case$log2fc <- as.numeric(case$log2fc)

fc_threshold <- sort(abs(case$log2fc), decreasing=T)[100]
case_filtered <- case[which(abs(case$log2fc) >= fc_threshold), ]

case_sorted <- case_filtered[order(case_filtered$log2fc), ]

case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene <- paste0( "***   ", case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene)

my_col <- rep("gray", dim(case_sorted)[2])
my_col[which(case_sorted$log2fc > 1)] <- "orange"
my_col[which(case_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(case_sorted$log2fc), na.rm=T)) + 3


pdf(paste0("barplot_", output_name, "_all_top100.pdf"))

barplot(case_sorted$log2fc, names.arg=case_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2fc", xlim=c(-my_xmax, my_xmax), cex.names=0.3, main=paste0(dim(case_filtered)[1], " proteins"))

abline(v=fc_threshold, lty=3)
abline(v=-fc_threshold, lty=3)
dev.off()

}


list_tsv <- system(command="ls *_preprocessed.tsv", intern=T)

cosmic <- fread("~/project/phrt/TpViz/docker-tp-viz/data/Cancer_Census_all_072018_COSMIC.csv")
list_prot <- fread("~/project/phrt/TpViz/docker-tp-viz/data/Melanoma_marker_Anja_IDs.csv", header=T)

for(i in 1:length(list_tsv)) {

filename <- gsub("_preprocessed.tsv", "", list_tsv[i])

message("Processing ", i, ": ", filename, "...")

raw <- fread(list_tsv[i], dec=",")

raw$log2fc <- as.numeric(raw$log2fc)

raw$Gene <- sapply( strsplit(raw$Gene, "\\|"), "[[", 1)

case <- copy(raw)

fc_threshold <- sort(abs(case$log2fc), decreasing=T)[100]
case_filtered <- case[which(abs(case$log2fc) >= fc_threshold), ]

case_sorted <- case_filtered[order(case_filtered$log2fc), ]

case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene <- paste0( "***   ", case_sorted[which(case_sorted$Gene %in% cosmic$Gene_Symbol), ]$Gene)

my_col <- rep("gray", dim(case_sorted)[2])
my_col[which(case_sorted$log2fc > 1)] <- "orange"
my_col[which(case_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(case_sorted$log2fc), na.rm=T)) + 3


pdf(paste0("barplot_", filename, "_top100.pdf"))

barplot(case_sorted$log2fc, names.arg=case_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2fc", xlim=c(-my_xmax, my_xmax), cex.names=0.3, main=paste0(dim(case_sorted)[1], " proteins"))

abline(v=fc_threshold, lty=3)
abline(v=-fc_threshold, lty=3)
dev.off()




markers <- raw[which(raw$Entry %in% list_prot$Uniprot), ]

markers$log2fc <- as.numeric(markers$log2fc)

markers_sorted <- markers[order(markers$log2fc), ]

my_col <- rep("gray", dim(markers_sorted)[1])
my_col[which(markers_sorted$log2fc > 1)] <- "orange"
my_col[which(markers_sorted$log2fc < -1)] <- "blue"

my_xmax <- round(max(abs(markers_sorted$log2fc), na.rm=T)) + 1


pdf(paste0("barplot_", filename, "_marker.pdf"))

barplot(markers_sorted$log2fc, names.arg=markers_sorted$Gene, horiz=T, las=1, col=my_col, xlab="log2FC", xlim=c(-my_xmax, my_xmax), cex.names=0.8, main=paste0(dim(markers_sorted)[1], " marker proteins"))
abline(v=1, lty=3)
abline(v=-1, lty=3)

dev.off()

message("Done.")

}
