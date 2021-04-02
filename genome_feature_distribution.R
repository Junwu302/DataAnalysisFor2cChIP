# library(GenomicRanges)
# library(GenomicAlignments)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(ChIPpeakAnno)
# bams = list.files(".",pattern = "*Rep.*rmdup.bam$")
# Reads_features = list()
# for(bam in bams){
#   print(paste0(c("Processing", bam, "..."),collapse = " "))
#   sampleA = granges(readGAlignmentPairs(bam))
#   aCR = assignChromosomeRegion(sampleA, nucleotideLevel=TRUE, 
#                                precedence=c("Promoters","fiveUTRs", "threeUTRs","Exons", "Introns"),
#                                TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
#                                proximal.promoter.cutoff=2000)
#   Reads_features[[gsub(".rmdup.bam","",bam)]] = aCR
# }
# save(Reads_features, file="Reads_features.RData")

library(reshape2)
library(ggplot2)
library(ggsci)
load("data/Reads_features.RData")
Reads_features = Reads_features[-c(2,4,13,16)]
df = do.call(rbind, lapply(Reads_features, function(x){
  x = x[[1]]
  x = 100*x/sum(x)
  return(x)
}))
df = melt(df)
df$Var1 = gsub("K562_H3K4me3_","",gsub("Rep","#",df$Var1))
df$Var2 = factor(df$Var2, levels = c("Intergenic.Region","immediateDownstream","threeUTRs",
                                     "Introns","Exons","fiveUTRs","Promoters"))
df$Var1 = factor(df$Var1, levels = c("10_#1","10_#3","40_#2","50_#1","100_#1","100_#2",
                                     "1000_#1","1000_#2","10-4_#1","10-4_#2","10-7_#1","10-7_#2"))
pdf("Reads_genome_feature.pdf",width = 10, height = 6)
ggplot(data = df, mapping = aes(x=Var1,y=value, fill=Var2)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_lancet() + ggtitle("K562 H3K4me3") + theme_bw()+ 
  scale_y_continuous(breaks=seq(0,100,20), expand = c(.005, 0),
                     labels=paste(seq(0,100,20),"%",sep='')) +
  guides(fill=guide_legend(reverse=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size=22, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=20, color = "black"),
        legend.position = "bottom")
dev.off()



