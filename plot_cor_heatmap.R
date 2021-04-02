library(pheatmap)
library(RColorBrewer)
cor_df = read.table("data/K562_H3K4me3_ReadCount_per_2k.tsv",sep='\t',header = T,stringsAsFactors = F)
colnames(cor_df)[-c(1:3)] = gsub(".rmdup.bam","",colnames(cor_df)[-c(1:3)])
cor_df = cor_df[,-c(1:4,9,12,15,18:21)]
colnames(cor_df) = gsub("(^[A-Za-z0-9]*_|cells|rep)","",colnames(cor_df))
cor_df = cor_df[,order(colnames(cor_df))]
mat = matrix(1, nrow = ncol(cor_df), ncol = ncol(cor_df))
for(i in 1:(ncol(cor_df)-1)){
  for(j in (i+1):ncol(cor_df)){
    df = cor_df[,c(i,j)]
    df = df[rowSums(df)>0,]
    mat[i,j] = mat[j,i] = cor(df[,1],df[,2])
  }
}
rownames(mat) = colnames(mat) = colnames(cor_df)
pdf("figures/K562_H3K4me3_Cor_heatmap.pdf")
breaksList = seq(0,1,0.1)
pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,display_numbers=T,number_color ="black")
dev.off()

cor_df = read.table("data/K562_H3K27ac_ReadCount_per_2k.tsv",sep='\t',header = T,stringsAsFactors = F)
colnames(cor_df)[-c(1:3)] = gsub(".rmdup.bam","",colnames(cor_df)[-c(1:3)])
cor_df = cor_df[,-c(1:4)]
colnames(cor_df) = gsub("(^[A-Za-z0-9]*_|cells|rep)","",colnames(cor_df))
cor_df = cor_df[,order(colnames(cor_df))]
mat = matrix(1, nrow = ncol(cor_df), ncol = ncol(cor_df))
for(i in 1:(ncol(cor_df)-1)){
  for(j in (i+1):ncol(cor_df)){
    df = cor_df[,c(i,j)]
    df = df[rowSums(df)>0,]
    mat[i,j] = mat[j,i] = cor(df[,1],df[,2])
  }
}
rownames(mat) = colnames(mat) = colnames(cor_df)
pdf("figures/K562_H3K27ac_Cor_heatmap.pdf")
breaksList = seq(0,1,0.1)
pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,display_numbers=T,number_color ="black")
dev.off()

cor_df = read.table("data/CUT_RUN/CUT_RUN_ReadCount_per_2k.tsv",sep='\t',header = T,stringsAsFactors = F)
colnames(cor_df)[-c(1:3)] = gsub(".rmdup.bam","",colnames(cor_df)[-c(1:3)])
cor_df = cor_df[,grep("H3",colnames(cor_df))]
colnames(cor_df) = gsub("(^[A-Za-z0-9]*_|cells|rep)","",colnames(cor_df))
cor_df = cor_df[,order(colnames(cor_df))]
mat = matrix(1, nrow = ncol(cor_df), ncol = ncol(cor_df))
for(i in 1:(ncol(cor_df)-1)){
  for(j in (i+1):ncol(cor_df)){
    df = cor_df[,c(i,j)]
    df = df[rowSums(df)>0,]
    mat[i,j] = mat[j,i] = cor(df[,1],df[,2])
  }
}
rownames(mat) = colnames(mat) = colnames(cor_df)
pdf("figures/CUT_RUN_Cor_heatmap.pdf")
breaksList = seq(0,1,0.1)
pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,display_numbers=T,number_color ="black")
dev.off()

cor_df = read.table("data/ChILSeq/ChILSeq_ReadCount_per_2k.tsv",sep='\t',header = T,stringsAsFactors = F)
colnames(cor_df)[-c(1:3)] = gsub(".rmdup.bam","",colnames(cor_df)[-c(1:3)])
cor_df = cor_df[,grep("H3",colnames(cor_df))]
colnames(cor_df) = gsub("(^[A-Za-z0-9]*_|cells|rep)","",colnames(cor_df))
cor_df = cor_df[,order(colnames(cor_df))]
mat = matrix(1, nrow = ncol(cor_df), ncol = ncol(cor_df))
for(i in 1:(ncol(cor_df)-1)){
  for(j in (i+1):ncol(cor_df)){
    df = cor_df[,c(i,j)]
    df = df[rowSums(df)>0,]
    mat[i,j] = mat[j,i] = cor(df[,1],df[,2])
  }
}
rownames(mat) = colnames(mat) = colnames(cor_df)
pdf("figures/ChILSeq_Cor_heatmap.pdf")
breaksList = seq(0,1,0.1)
pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,display_numbers=T,number_color ="black")
dev.off()



###########################

pairs = combn(colnames(cor_df), 2)
for(i in 1:ncol(pairs)){
  S1 = pairs[1,i]
  S2 = pairs[2,i]
  df = cor_df[,c(S1, S2)]
  R = round(cor(df[,1], df[,2]),3)
  df = 1e6*df/colSums(df)/2
  df = df[apply(df, 1, min)>0,]
  df = log10(df)
  colnames(df) = c("x","y")
  pdf(paste0(c("figures/Cor_Scatter_Fig/",S1,"_",S2,"_Cor_Scatter.pdf"),collapse = ""))
  smoothScatter(df$x,df$y,xlab = paste(S1,"log10(RPKM)"), ylab = paste(S2,"log10(RPKM)"),
                colramp =colorRampPalette(c("white","blue", "green","orange", "red")))
  text(x= min(df$x)+(max(df$x)-min(df$x))/10,
       y=max(df$y)-(max(df$y)-min(df$y))/10,labels = paste("R = ", R,sep=""),col="darkred")
  dev.off()
}

### TSS heatmap
TSS_mat = read.table("data/K562_H3K4me3_TSS_matrix",skip = 1,sep='\t',stringsAsFactors = F)
TSS_mat_new = TSS_mat[,1:6]
TSS_mat = TSS_mat[,-c(1:6)]
for(i in 1:12){
  df = TSS_mat[,(120*(i-1)+1):(120*i)]
  df = round(df/sum(df)*1e6,4)
  TSS_mat_new = cbind(TSS_mat_new, df)
}
write.table(TSS_mat_new, file= "K562_H3K4me3_TSS_matrix_new",sep='\t',row.names = F, col.names = F, quote = F)

### ENCFF752ALB heatmap
ENCFF752ALB_mat = read.table("data/K562_H3K4me3_ENCFF752ALB_matrix",skip = 1,sep='\t',stringsAsFactors = F)
ENCFF752ALB_mat_new = ENCFF752ALB_mat[,1:6]
ENCFF752ALB_mat = ENCFF752ALB_mat[,-c(1:6)]
for(i in 1:12){
  df = ENCFF752ALB_mat[,(120*(i-1)+1):(120*i)]
  df = round(df/sum(df)*1e6,4)
  ENCFF752ALB_mat_new = cbind(ENCFF752ALB_mat_new, df)
}
write.table(ENCFF752ALB_mat_new, file= "K562_H3K4me3_ENCFF752ALB_matrix_new",sep='\t',row.names = F, col.names = F, quote = F)
