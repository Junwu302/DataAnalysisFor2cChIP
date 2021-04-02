library(GenomicRanges)
load_MACS2_peak <- function(file){
  require(GenomicRanges, quietly = T)
  peak = read.table(file, sep='\t',header = F, stringsAsFactors = F)
  peak = peak[,c(1:5,9)]
  colnames(peak) = c("chrom","start","end","name","summit","qval")
  peak = peak[peak$chrom %in% paste("chr",c(1:22,"X","Y"),sep=''),]
  peak$qval = 10^(-peak$qval)
  peak$summit = peak$start + peak$summit
  peak = GRanges(seqnames = Rle(peak$chrom),ranges = IRanges(start = peak$start, end = peak$end),
                 name = peak$name, summit = peak$summit, qval = peak$qval)
  return(peak)
}
get.ROC <- function(queryPeak, refPeak, anchor){
  require(plyr, quietly = T)
  require(pROC, quietly = T)
  anchor$name = paste("Feature",1:length(anchor),sep=':')
  gold_standard = overlapsAny(anchor, refPeak)
  df = data.frame(findOverlaps(anchor, queryPeak), stringsAsFactors = F)
  df$queryHits = anchor$name[df$queryHits]
  df$subjectHits = queryPeak$qval[df$subjectHits]
  df = ddply(df, .variables = "queryHits",.fun = function(x){min(x$subjectHits)})
  pred = rep(1, length(anchor))
  names(pred) = anchor$name
  pred[df$queryHits] = df$V1
  res = roc(response=factor(gold_standard), predictor = pred, percent=T)
  return(res)
}

load("data/hg19_promoter.RData") #(+2000, -500)
hg19_promoter = promoters(hg19_promoter, upstream = 0, downstream = 4000)
K562_H3K4me3_RefPeak = load_MACS2_peak("data/Peak_for_ROC/K562_H3K4me3_ENCFF752ALB.bed")
K562_H3K4me3_RefPeak = K562_H3K4me3_RefPeak[K562_H3K4me3_RefPeak$qval < 0.01,]
K562_H3K27ac_RefPeak = load_MACS2_peak("data/Peak_for_ROC/K562_H3K27ac_ENCFF044JNJ.bed")
K562_H3K27ac_RefPeak = K562_H3K27ac_RefPeak[K562_H3K27ac_RefPeak$qval < 0.01,]

## ROC analysis for K562 H3K4me3
peak_files = list.files("data/Peak_for_ROC/",pattern = "K562_H3K4me3.*_roc_peaks.*Peak$",full.names = T)
K562_H3K4me3_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(K562_H3K4me3_peaks) = gsub("_roc_peaks.*","",basename(peak_files))

K562_H3K4me3_public_ROC = lapply(K562_H3K4me3_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},K562_H3K4me3_RefPeak, hg19_promoter)
for(i in 1:length(K562_H3K4me3_public_ROC)){
  print(paste(names(K562_H3K4me3_public_ROC)[i], round(auc(K562_H3K4me3_public_ROC[[i]]),3),sep = ':'))
}
save(K562_H3K4me3_public_ROC, file="data/K562_H3K4me3_public_ROC.RData")
K562_H3K4me3_public_ROC = lapply(K562_H3K4me3_peaks[-c(4:6)], function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},K562_H3K4me3_peaks[[6]], hg19_promoter)
for(i in 1:length(K562_H3K4me3_public_ROC)){
  print(paste(names(K562_H3K4me3_public_ROC)[i], round(auc(K562_H3K4me3_public_ROC[[i]]),3),sep = ':'))
}
save(K562_H3K4me3_public_ROC, file="data/K562_H3K4me3_public_ROC.RData")

## ROC analysis for K562 H3K27ac
peak_files = list.files("data/Peak_for_ROC/",pattern = "K562_H3K27ac.*_roc_peaks.*Peak$",full.names = T)
K562_H3K27ac_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(K562_H3K27ac_peaks) = gsub("_roc_peaks.*","",basename(peak_files))

K562_H3K27ac_public_ROC = lapply(K562_H3K27ac_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},K562_H3K27ac_RefPeak, hg19_promoter)
for(i in 1:length(K562_H3K27ac_public_ROC)){
  print(paste(names(K562_H3K27ac_public_ROC)[i], round(auc(K562_H3K27ac_public_ROC[[i]]),3),sep = ':'))
}
save(K562_H3K27ac_public_ROC, file="data/K562_H3K27ac_public_ROC.RData")


## ROC analysis for CUT&RUN
load("data/mm10_promoter.RData") #(+2000, -500)
mm10_promoter = promoters(mm10_promoter, upstream = 0, downstream = 4000)
ESC_H3K4me3_RefPeak = load_MACS2_peak("data/CUT_RUN/ENCFF204SXD_H3K4me3_peaks.narrowPeak")
ESC_H3K4me3_RefPeak = ESC_H3K4me3_RefPeak[ESC_H3K4me3_RefPeak$qval < 0.01,]
peak_files = list.files("data/CUT_RUN/",pattern = "*H3K4me3.*_roc_peaks.*Peak$",full.names = T)
ESC_H3K4me3_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(ESC_H3K4me3_peaks) = gsub("_roc_peaks.*","",basename(peak_files))
CUT_RUN_ESC_H3K4me3_public_ROC = lapply(ESC_H3K4me3_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},ESC_H3K4me3_RefPeak, mm10_promoter)
for(i in 1:length(CUT_RUN_ESC_H3K4me3_public_ROC)){
  print(paste(names(CUT_RUN_ESC_H3K4me3_public_ROC)[i], round(auc(CUT_RUN_ESC_H3K4me3_public_ROC[[i]]),3),sep = ':'))
}
save(CUT_RUN_ESC_H3K4me3_public_ROC, file="data/CUT_RUN_ESC_H3K4me3_public_ROC.RData")

ESC_H3K27ac_RefPeak = load_MACS2_peak("data/CUT_RUN/ENCFF587LEB_H3K27ac_peaks.narrowPeak")
ESC_H3K27ac_RefPeak = ESC_H3K27ac_RefPeak[ESC_H3K27ac_RefPeak$qval < 0.01,]
peak_files = list.files("data/CUT_RUN/",pattern = "*H3K27ac.*_roc_peaks.*Peak$",full.names = T)
ESC_H3K27ac_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(ESC_H3K27ac_peaks) = gsub("_roc_peaks.*","",basename(peak_files))
CUT_RUN_ESC_H3K27ac_public_ROC = lapply(ESC_H3K27ac_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},ESC_H3K27ac_RefPeak, mm10_promoter)
for(i in 1:length(CUT_RUN_ESC_H3K27ac_public_ROC)){
  print(paste(names(CUT_RUN_ESC_H3K27ac_public_ROC)[i], round(auc(CUT_RUN_ESC_H3K27ac_public_ROC[[i]]),3),sep = ':'))
}
save(CUT_RUN_ESC_H3K27ac_public_ROC, file="data/CUT_RUN_ESC_H3K27ac_public_ROC.RData")


## ROC analysis for ChILSeq
load("data/mm10_promoter.RData") #(+2000, -500)
mm10_promoter = promoters(mm10_promoter, upstream = 0, downstream = 4000)
C2C12_H3K4me3_RefPeak = load_MACS2_peak("data/ChILSeq/C2C12_H3K4me3_million_peaks.narrowPeak")
C2C12_H3K4me3_RefPeak = C2C12_H3K4me3_RefPeak[C2C12_H3K4me3_RefPeak$qval < 0.01,]
peak_files = list.files("data/ChILSeq/",pattern = "*H3K4me3.*_roc_peaks.*Peak$",full.names = T)
C2C12_H3K4me3_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(C2C12_H3K4me3_peaks) = gsub("_roc_peaks.*","",basename(peak_files))
ChILSeq_C2C12_H3K4me3_public_ROC = lapply(C2C12_H3K4me3_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},C2C12_H3K4me3_RefPeak, mm10_promoter)
for(i in 1:length(ChILSeq_C2C12_H3K4me3_public_ROC)){
  print(paste(names(ChILSeq_C2C12_H3K4me3_public_ROC)[i], round(auc(ChILSeq_C2C12_H3K4me3_public_ROC[[i]]),3),sep = ':'))
}
save(ChILSeq_C2C12_H3K4me3_public_ROC, file="data/ChILSeq_C2C12_H3K4me3_public_ROC.RData")


C2C12_H3K27ac_RefPeak = load_MACS2_peak("data/ChILSeq/C2C12_H3K27ac_million_peaks.narrowPeak")
C2C12_H3K27ac_RefPeak = C2C12_H3K27ac_RefPeak[C2C12_H3K27ac_RefPeak$qval < 0.01,]
peak_files = list.files("data/ChILSeq/",pattern = "*H3K27ac.*_roc_peaks.*Peak",full.names = T)
C2C12_H3K27ac_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(C2C12_H3K27ac_peaks) = gsub("_roc_peaks.*","",basename(peak_files))
ChILSeq_C2C12_H3K27ac_public_ROC = lapply(C2C12_H3K27ac_peaks, function(peak, bulk, anchor){
  return(get.ROC(peak, bulk, anchor))
},C2C12_H3K27ac_RefPeak, mm10_promoter)
for(i in 1:length(ChILSeq_C2C12_H3K27ac_public_ROC)){
  print(paste(names(ChILSeq_C2C12_H3K27ac_public_ROC)[i], round(auc(ChILSeq_C2C12_H3K27ac_public_ROC[[i]]),3),sep = ':'))
}
save(ChILSeq_C2C12_H3K27ac_public_ROC, file="data/ChILSeq_C2C12_H3K27ac_public_ROC.RData")




### Plot
library(ggplot2)
library(RColorBrewer)
roc_files = list.files("data/",pattern = "*ROC.RData$",full.names = T)
for(fl in roc_files){
  load(fl)
}
###================== H3K4me3 ROC Plot ===============================
## H3k4me3 1K cells
df1 = list(AHCX_H3K4me3_1K_rep1 = K562_H3K4me3_public_ROC[[12]],
           AHCX_H3K4me3_1K_rep2 = K562_H3K4me3_public_ROC[[13]],
           ChIL_H3K4me3_1K_rep1 = ChILSeq_C2C12_H3K4me3_public_ROC[[1]],
           ChIL_H3K4me3_1K_rep2 = ChILSeq_C2C12_H3K4me3_public_ROC[[2]],
           CUT_H3K4me3_500_rep1 = CUT_RUN_ESC_H3K4me3_public_ROC[[1]],
           CUT_H3K4me3_500_rep2 = CUT_RUN_ESC_H3K4me3_public_ROC[[4]])
H3K4me3_1K_auc = unlist(lapply(df1, function(x){auc(x)}))
df1 = do.call(rbind, lapply(df1, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df1$Sample = gsub("\\.[0-9]*$","",rownames(df1))
pdf("figures/H3K4me3_1K_ROC.pdf", width = 8, height = 4)
ggplot(data = df1, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K4me3_1K_auc)), round(H3K4me3_1K_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2")+ theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()
## H3k4me3 100 cells
df2 = list(AHCX_H3K4me3_100_rep1 = K562_H3K4me3_public_ROC[[9]],
           AHCX_H3K4me3_100_rep2 = K562_H3K4me3_public_ROC[[10]],
           ChIL_H3K4me3_100_rep1 = ChILSeq_C2C12_H3K4me3_public_ROC[[3]],
           ChIL_H3K4me3_100_rep2 = ChILSeq_C2C12_H3K4me3_public_ROC[[4]])
H3K4me3_100_auc = unlist(lapply(df2, function(x){auc(x)}))
df2 = do.call(rbind, lapply(df2, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df2$Sample = gsub("\\.[0-9]*$","",rownames(df2))
pdf("figures/H3K4me3_100_ROC.pdf", width = 8, height = 4)
ggplot(data = df2, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K4me3_100_auc)), round(H3K4me3_100_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2")+ theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
      axis.title = element_text(colour = "black",size = 16),
      legend.title=element_blank(),
      legend.text = element_text(colour = "black",size=14))
dev.off()
## H3k4me3 50 cells
df3 = list(AHCX_H3K4me3_50_rep1 = K562_H3K4me3_public_ROC[[18]],
           AHCX_H3K4me3_50_rep2 = K562_H3K4me3_public_ROC[[19]],
           CUT_H3K4me3_50_rep1 = CUT_RUN_ESC_H3K4me3_public_ROC[[2]],
           CUT_H3K4me3_50_rep2 = CUT_RUN_ESC_H3K4me3_public_ROC[[5]])
H3K4me3_50_auc = unlist(lapply(df3, function(x){auc(x)}))
df3 = do.call(rbind, lapply(df3, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df3$Sample = gsub("\\.[0-9]*$","",rownames(df3))
pdf("figures/H3K4me3_50_ROC.pdf", width = 8, height = 4)
ggplot(data = df3, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K4me3_50_auc)), round(H3K4me3_50_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()
## H3k4me3 10 cells
df4 = list(AHCX_H3K4me3_10_rep1 = K562_H3K4me3_public_ROC[[4]],
           AHCX_H3K4me3_10_rep2 = K562_H3K4me3_public_ROC[[5]],
           CUT_H3K4me3_10_rep1 = CUT_RUN_ESC_H3K4me3_public_ROC[[3]],
           CUT_H3K4me3_10_rep2 = CUT_RUN_ESC_H3K4me3_public_ROC[[6]])
H3K4me3_10_auc = unlist(lapply(df4, function(x){auc(x)}))
df4 = do.call(rbind, lapply(df4, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df4$Sample = gsub("\\.[0-9]*$","",rownames(df4))
pdf("figures/H3K4me3_10_ROC.pdf", width = 8, height = 4)
ggplot(data = df4, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K4me3_10_auc)), round(H3K4me3_10_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()
###================== H3K27ac ROC Plot=====================
df5 = list(AHCX_H3K27ac_1K_rep1 = K562_H3K27ac_public_ROC[[5]],
           AHCX_H3K27ac_1K_rep2 = K562_H3K27ac_public_ROC[[6]],
           ChIL_H3K27ac_1K_rep1 = ChILSeq_C2C12_H3K27ac_public_ROC[[1]],
           ChIL_H3K27ac_1K_rep2 = ChILSeq_C2C12_H3K27ac_public_ROC[[2]])
H3K27ac_1K_auc = unlist(lapply(df5, function(x){auc(x)}))
df5 = do.call(rbind, lapply(df5, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df5$Sample = gsub("\\.[0-9]*$","",rownames(df5))
pdf("figures/H3K27ac_1K_ROC.pdf", width = 8, height = 4)
ggplot(data = df5, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K27ac_1K_auc)), round(H3K27ac_1K_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()
  
df6 = list(AHCX_H3K27ac_100_rep1 = K562_H3K27ac_public_ROC[[3]],
           AHCX_H3K27ac_100_rep2 = K562_H3K27ac_public_ROC[[4]],
           ChIL_H3K27ac_100_rep1 = ChILSeq_C2C12_H3K27ac_public_ROC[[3]],
           ChIL_H3K27ac_100_rep2 = ChILSeq_C2C12_H3K27ac_public_ROC[[4]])
H3K27ac_100_auc = unlist(lapply(df6, function(x){auc(x)}))
df6 = do.call(rbind, lapply(df6, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df6$Sample = gsub("\\.[0-9]*$","",rownames(df6))
pdf("figures/H3K27ac_100_ROC.pdf", width = 8, height = 4)
ggplot(data = df6, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K27ac_100_auc)), round(H3K27ac_100_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()

df7 = list(AHCX_H3K27ac_50_rep1 = K562_H3K27ac_public_ROC[[7]],
           AHCX_H3K27ac_50_rep2 = K562_H3K27ac_public_ROC[[8]],
           CUT_H3K27ac_50_rep1 = CUT_RUN_ESC_H3K27ac_public_ROC[[1]],
           CUT_H3K27ac_50_rep2 = CUT_RUN_ESC_H3K27ac_public_ROC[[2]])
H3K27ac_50_auc = unlist(lapply(df7, function(x){auc(x)}))
df7 = do.call(rbind, lapply(df7, function(x){
  data.frame(sensitivity=x$sensitivities, specificity=x$specificities, stringsAsFactors = F)
}))
df7$Sample = gsub("\\.[0-9]*$","",rownames(df7))
pdf("figures/H3K27ac_50_ROC.pdf", width = 8, height = 4)
ggplot(data = df7, mapping = aes(x=specificity, y =sensitivity)) +
  geom_line(mapping = aes(color=Sample),size=1) + xlim(100,0) + ylim(0,100) +
  ggtitle(paste0(paste(gsub("_.*","",names(H3K27ac_50_auc)), round(H3K27ac_50_auc,3)),collapse = ";")) + 
  scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(axis.text = element_text(colour = "black",size=14),
        axis.title = element_text(colour = "black",size = 16),
        legend.title=element_blank(),
        legend.text = element_text(colour = "black",size=14))
dev.off()

df8 = list(AHCX_H3K27ac_10_rep1 = K562_H3K27ac_public_ROC[[1]],
           AHCX_H3K27ac_10_rep2 = K562_H3K27ac_public_ROC[[2]])

