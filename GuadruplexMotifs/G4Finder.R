dist <- data.frame(read.table("/Users/zeidh/Downloads/nullDistribution.txt"))
dist2 <- data.frame(read.table("/Users/zeidh/Desktop/Enrichment/nullDistribution.txt"))
dist3 <- data.frame(read.table("/Users/zeidh/Desktop/Enrichment/exonDistribution.txt"))
library(tidyverse)

str(dist2)
ggplot(dist2) + geom_density(aes(V1),
							fill = "darkgrey") +
	geom_vline(aes(xintercept = dist2[1,]), 
			   linetype = "dashed", size = 0.6,
			   color = "#FC4E07") +
	theme_classic() + ggsave("promoterEnrichment.jpg")


library(pqsfinder)

if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
BiocManager::install("pqsfinder")

seq <- DNAString("TTTTGGGAGGGTGGGAGGGA")
pqsfinder(seq, min_score = 53, run_min_len=3,run_max_len=3,loop_min_len = 0, loop_max_len = 7)

library(pqsfinder)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
BiocManager::install("rtracklayer")
library(rtracklayer)
library(ggplot2)
BiocManager::install("Gviz")
library(Gviz)

hg <- system.file("extdata","Hsapiens.fa",package="BSgenome")


gnm <- "hg38"

prefix="chr"
suffix = c(1:21,"X","Y")
chromosomes <- paste(prefix,suffix,sep="")

for (chr in chromosomes){
	seq <- Hsapiens[[chr]]
	pqs <- pqsfinder(seq, min_score = 20, run_min_len=3,run_max_len=3,loop_min_len = 0, loop_max_len = 7)
	name <- paste0(chr,"pqs")
	assign(name,pqs)
}


suffix = "pqs"
chromosomeObjects <- paste(chromosomes,suffix,sep="")
chromosomes
bed<-data.frame()
count=0
for (obj in chromosomes){
	obj=chromosomes[1]
	name <- paste0(obj,"pqs")
	object <- get(name)
	bed <- data.frame(chromosome=obj,start=object@ranges@start,stop=(object@ranges@start+object@ranges@width),width=object@ranges@width,score=object@elementMetadata@listData$score, strand=object@elementMetadata@listData$strand)
	write.table(bed,"G4s.bed",append=T,sep="\t",row.names = F,col.names = F,quote = F)
	count = count + nrow(bed)
}
print(count)



gains <- toGRanges(data.frame(chr=c("chr1", "chr5", "chr17", "chr22"), start=c(1, 1000000, 4000000, 1),
							  end=c(5000000, 3200000, 80000000, 1200000)))
losses <- toGRanges(data.frame(chr=c("chr3", "chr9", "chr17"), start=c(80000000, 20000000, 1),
							   end=c(170000000, 30000000, 25000000)))
kp <- plotKaryotype(genome="hg19")
kpPlotRegions(kp, gains, col="#FFAACC")
kpPlotRegions(kp, losses, col="#CCFFAA")

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chr1", "chr2", "chr3"))
