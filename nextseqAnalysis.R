nextseq <- data.frame(read.table("/Users/zeidh/Desktop/Feb25NextSeq.txt",header=T,fill=T,sep="\t"))
library(tidyverse)
nextseq <- nextseq %>% select(-c(X,X.1,X.2)) 
nextseq$library = nextseq$Library
nextseq <- nextseq[1:72,]
str(nextseq)

library(stringr)
library(dplyr)
library(ggpubr)

nextseq <- nextseq %>%
	separate(library, c("name", "gene","UV","hoescht","cluster","well","S"), "_")

nextseq <- nextseq[!is.na(nextseq$SCE_down),]


nextseq$gene <- as.factor(nextseq$gene)
nextseq$SCE_down <- as.numeric(nextseq$SCE_down)
nextseq$SCE <- as.numeric(nextseq$SCE_total)
nextseq$SCE <- as.numeric(nextseq$SCE_big)
nextseq$SCE <- as.numeric(nextseq$SCE_small)
nextseq$Coverage <- as.numeric(sub("%", "", nextseq$Coverage))

nextseq$SCE_norm_ploidy = nextseq$SCE_down/nextseq$ploidy
nextseq$SCE_norm_reads = (nextseq$SCE_norm_ploidy/nextseq$Reads_aligned_postfiltering)



(plot <- ggplot(nextseq) + 
	geom_point(aes(Reads_aligned_postfiltering,SCE_small, color=gene)) + 
	geom_smooth(se=F,aes(Reads_aligned_postfiltering,SCE_small, color=gene)) +
	theme_classic() + 
	labs(x="Reads Aligned Post Filtering",y="Sister Chromatid Exchange Count",title="SCE Count vs Reads Count") +
	theme(legend.title = element_blank())) # +
	#ggsave("SCEvsReadCount(outliersRemoved).png"))

(plot2 <- ggplot(nextseq) + 
		geom_point(aes(SCE,Coverage, color=gene)) + 
		geom_smooth(se=F,aes(SCE,Coverage, color=gene)) +
		theme_classic() + 
		labs(x="Sister Chromatid Exchange Count",y="Coverage",title="SCE Count vs Coverage") +
		theme(legend.title = element_blank())  +
		ggsave("SCEvsCoverage.png"))

(plot3 <- ggplot(filteredTwice) + 
		geom_point(aes(SCE,Sequencing_for_5_percent_coverage, color=gene)) + 
		geom_smooth(method = loess, method.args = list(family = "symmetric"),se=F,aes(SCE,Sequencing_for_5_percent_coverage)) +
		theme_classic() + 
		labs(x="Sister Chromatid Exchange Count",y="Sequencing_for_5_percent_coverage",title="SCE Count vs Complexity") +
		theme(legend.title = element_blank())  +
		ggsave("SCEvsComplexity(outliersRemoved).png"))

(plot4 <- ggplot(nextseq) +
		geom_violin(aes(gene, SCE_total,color=gene,fill=gene))+
		geom_jitter(aes(gene,SCE_total,color=gene))+
		scale_fill_manual(values=c("#7a281b", "#195e23", "#3943ad"))+
		theme_classic() + 
		labs(y="Sister Chromatid Exchange Count",x="Gene Knockout",title="Sister Chromatid Exchange Events") +
		theme(legend.title = element_blank()))#  +
	#	ggsave("SCEvsGene.png"))

(plot5<- ggplot(nextseq) +
		geom_boxplot(aes(gene, norm_sce,color=gene,fill=gene))+
		geom_jitter(aes(gene,norm_sce,color=gene))+
		scale_fill_manual(values=c("#7a281b", "#195e23", "#3943ad"))+
		theme_classic() + 
		labs(y="Sister Chromatid Exchange Count",x="Gene Knockout",title="Sister Chromatid Exchange Events") +
		theme(legend.title = element_blank()))

nextseq$norm_sce  <- ((nextseq$SCE_big/nextseq$Reads_aligned_postfiltering)/nextseq$ploidy)


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

ggboxplot(nextseq,"gene", "SCE_norm_reads",
		  color="gene",
		  palette = c("#0BAE7C", "#D88907", "#AE0B3D"),
		  #fill= c("#85D6BE","#E7B564","#CC6786"),
		  add = "jitter",add.params = list(color="gene")) +
	stat_compare_means(method = "anova", label.y = 3e-5) +      # Add global p-value
	stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")# +
#ggsave("initialSCE_count.png")# Pairwise comparison against reference     # Add global p-value



compare_means(SCE_norm_ploidy ~ gene,  data = nextseq, method = "anova")
res.aov <- aov(SCE_norm_ploidy~ gene, data = nextseq)
summary(res.aov)
TukeyHSD(res.aov)

nextseq[order(-nextseq$Reads_aligned_postfiltering),]
count(nextseq, gene)

nextseq[order(nextseq$Reads_aligned_prefiltering),]
