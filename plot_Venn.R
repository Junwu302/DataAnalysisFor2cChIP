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
get_venn <- function(peaks, anchor, gold_peak){
  A = overlapsAny(anchor, gold_peak)
  venn_res = data.frame(do.call(rbind,lapply(peaks, function(peak, anchor, A){
    B = overlapsAny(anchor, peak)
    coverage = sum(A & B)/sum(A)
    preciesion = sum(A & B)/sum(B)
    return(c(preciesion = preciesion, coverage = coverage))
  }, anchor, A)))
  return(venn_res)
}

anchor_list = list(promoter=hg19_promoter, bin1k = hg19_1Kb_bins,bin4k = hg19_4Kb_bins)

load("data/hg19_promoter.RData") #(+2000, -500)
hg19_promoter = promoters(hg19_promoter, upstream = 0, downstream = 4000)
K562_H3K4me3_RefPeak = load_MACS2_peak("data/Peak_for_ROC/K562_H3K4me3_ENCFF752ALB.bed")
K562_H3K4me3_RefPeak = K562_H3K4me3_RefPeak[K562_H3K4me3_RefPeak$qval < 0.01,]
K562_H3K27ac_RefPeak = load_MACS2_peak("data/Peak_for_ROC/K562_H3K27ac_ENCFF044JNJ.bed")
K562_H3K27ac_RefPeak = K562_H3K27ac_RefPeak[K562_H3K27ac_RefPeak$qval < 0.01,]

peak_files = list.files("data/K562_peaks/",pattern = "K562_H3K4me3.*_peaks.*Peak$",full.names = T)
K562_H3K4me3_peaks = lapply(peak_files, function(fl){
  load_MACS2_peak(fl)
})
names(K562_H3K4me3_peaks) = gsub("_peaks.*","",basename(peak_files))
K562_H3K4me3_Venn = get_venn(K562_H3K4me3_peaks, hg19_promoter, K562_H3K4me3_RefPeak)




## H3k4me3
K562_H3K4me3_RefPeak = load_MACS2_peak("data/Peak_for_ROC/K562_H3K4me3_ENCFF752ALB.bed")
### 50 cell
fls = list.files("./H3K4me3_peak_PR/K562_H3K4me3_50_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K4me3_peak_PR/K562_H3K4me3_50_roc/",fl,sep=''))
})
K562_H3K4me3_venn_50cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K4me3_venn_50cell, file = "K562_H3K4me3_venn_50cell.RData")
### 100 cell
fls = list.files("./H3K4me3_peak_PR/K562_H3K4me3_100_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K4me3_peak_PR/K562_H3K4me3_100_roc/",fl,sep=''))
})
K562_H3K4me3_venn_100cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K4me3_venn_100cell, file = "K562_H3K4me3_venn_100cell.RData")
### 1000 cell
fls = list.files("./H3K4me3_peak_PR/K562_H3K4me3_1000_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K4me3_peak_PR/K562_H3K4me3_1000_roc/",fl,sep=''))
})
K562_H3K4me3_venn_1000cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K4me3_venn_1000cell, file = "K562_H3K4me3_venn_1000cell.RData")


## H3K27ac
## H3K27ac
gold_peak = load_peak("SRR6348898_H3K27ac_peaks.narrowPeak")
### 10 cell
fls = list.files("./H3K27ac_peak_PR/K562_H3K27ac_10_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K27ac_peak_PR/K562_H3K27ac_10_roc/",fl,sep=''))
})
K562_H3K27ac_venn_10cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K27ac_venn_10cell, file = "K562_H3K27ac_venn_10cell.Rata")
### 50 cell
fls = list.files("./H3K27ac_peak_PR/K562_H3K27ac_50_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K27ac_peak_PR/K562_H3K27ac_50_roc/",fl,sep=''))
})
K562_H3K27ac_venn_50cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K27ac_venn_50cell, file = "K562_H3K27ac_venn_50cell.Rata")
### 100 cell
fls = list.files("./H3K27ac_peak_PR/K562_H3K27ac_100_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K27ac_peak_PR/K562_H3K27ac_100_roc/",fl,sep=''))
})
K562_H3K27ac_venn_100cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K27ac_venn_100cell, file = "K562_H3K27ac_venn_100cell.Rata")
### 1000 cell
fls = list.files("./H3K27ac_peak_PR/K562_H3K27ac_1000_roc/",pattern = ".broadPeak$")
names(fls) = gsub("_peaks.broadPeak","", fls)
peaks = lapply(fls, function(fl){
  load_peak(paste("./H3K27ac_peak_PR/K562_H3K27ac_1000_roc/",fl,sep=''))
})
K562_H3K27ac_venn_1000cell = get_venn(peaks = peaks, anchor_list = anchor_list)
save(K562_H3K27ac_venn_1000cell, file = "K562_H3K27ac_venn_1000cell.Rata")




## ---------------------- get venn results ---------------------- ##

load("K562_H3K4me3_venn_100cell.RData")

