### Analyzing SCE data from Feb 25 NextSeq run

#required packages
library(tidyverse)
library(stringr)
library(dplyr)
library(ggpubr)

#load data
nextseq <- data.frame(read.table("/Users/zeidh/Desktop/Feb\ 25\ NextSeq/Feb25NextSeq.txt",header=T,fill=T,sep="\t"))


#clean up data
nextseq$library = nextseq$Library
nextseq <- nextseq %>%
	separate(library, c("name", "gene","UV","hoescht","cluster","well","S"), "_")
nextseq$gene <- as.factor(nextseq$gene)
nextseq$SCE_down <- as.numeric(nextseq$SCE_down)
nextseq$SCE <- as.numeric(nextseq$SCE_total)
nextseq$Coverage <- as.numeric(sub("%", "", nextseq$Coverage))
#nextseq <- nextseq %>% select(-c(X,X.1,X.2)) 
#nextseq <- nextseq[1:72,]
nextseq$SCE_norm_ploidy = nextseq$SCE_down/nextseq$ploidy
nextseq$SCE_norm_reads = (nextseq$SCE_norm_ploidy/nextseq$Reads_aligned_postfiltering)
nextseq$norm_sce  <- ((nextseq$SCE_big/nextseq$Reads_aligned_postfiltering)/nextseq$ploidy)
count(nextseq, gene)


#filter for good libraries
nextseq <- nextseq[!is.na(nextseq$SCE_down),]


#determine significance
res.aov <- aov(SCE_norm_ploidy~ gene, data = nextseq)
summary(res.aov)
TukeyHSD(res.aov)



#SCE count vs read count 
plot <- function(data,x,y,color){
	ggplot(data) + 
	geom_point(aes(x,y, color=color)) + 
	geom_smooth(se=F,aes(x,y,color=color)) +
	theme_classic() + 
	theme(legend.title = element_blank()) 
}


ggboxplot(nextseq, x = "gene", y = "norm_sce",
		  color = "gene",
		  palette = c("#0BAE7C", "#D88907", "#AE0B3D"),
		  fill= c("#85D6BE","#E7B564","#CC6786"),
		  add = "jitter") +
	stat_compare_means(method = "anova", label.y = 5e-05) +      # Add global p-value
	stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")  +
	ggsave("normalized_sce.png")# Pairwise comparison against reference     # Add global p-value

ggboxplot(nextseq,"gene", "SCE_norm_ploidy",
		color="gene",
		palette = c("#0BAE7C", "#D88907", "#AE0B3D"),
		#fill= c("#85D6BE","#E7B564","#CC6786"),
		add = "jitter",add.params = list(color="gene")) +
		stat_compare_means(method = "anova", label.y = 9
						) +      # Add global p-value
		stat_compare_means(label = "p.signif", method = "t.test",ref.group = "WT")# +
		ggsave("downSamp300k_ploidyNorm.png")# Pairwise comparison against reference     # Add global p-value




