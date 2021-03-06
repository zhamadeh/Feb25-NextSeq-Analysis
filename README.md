# Feb25 NextSeq sequencing run
On February 25th, we sequenced 24 libraries from three different samples: one *BLM-/-* line, one *RECQL5 -/-* line and one WT line with the intent of mapping sister chromatid exchange (SCE) events to the genome and compare models.

## Analysis 
Because of different sequencing depths among libraries, the following analysis takes a downsampling approach to normalize the level of reads to a constant level (350k reads in this case, but can be higher). Libraries were also normalized by ploidy (the WT line is haploid while the other two are diploid.

## Results

#### SCE event counts normalized by ploidy
![](https://github.com/zhamadeh/Feb25-NextSeq-Analysis/blob/master/Plots/unfilteredCounts.png)

#### SCE detection rate per 350k reads normalized by ploidy
![](https://github.com/zhamadeh/Feb25-NextSeq-Analysis/blob/master/Plots/downSamp300k_ploidyNorm.png)


