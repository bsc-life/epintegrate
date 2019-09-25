# LIBRARIES
library("dplyr")
library("ggplot2")
library("gridExtra")
library("cowplot")

colors <- c("#00AFBB", "#E7B800", "#FC4E07")
colors2 <- c("#FFDB6D","#D16103","#C3D7A4")

###
CCscore <- data.frame(all.scores[[1]][[1]])
colnames(CCscore) <- c("score")
plot.CC <- ggplot(CCscore, aes(rownames(CCscore),score, group=1)) + geom_line(color=colors[1]) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab("Chromatin contacts score")

###
# Histones
v.histones.chr1 <- read.delim("importance_vector_histones.csv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

bins.100kb <- data.frame(stringsAsFactors = FALSE)
for (i in seq(1,nrow(v.histones.chr1),5)){
  new.bin.val <- mean(v.histones.chr1[i:i+5,])
  bins.100kb <- rbind(bins.100kb,new.bin.val)
}
colnames(bins.100kb) <- c("score")

plot.hist <- ggplot(bins.100kb, aes(rownames(bins.100kb),score,group=1)) + geom_line(color=colors[3]) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab("Histone marks score")

###
# Add SNPs
snps.data <- read.delim("~/Documents/Scripts/Analysis/New/variants.data.csv", stringsAsFactors = FALSE)
snps.chr1 <- select(filter(snps.data, CHR_ID == 1), CHR_ID, CHR_POS)

get.bin <- function(x){
  return(as.character(floor(x/100000)))
}
snps.chr1$BIN <- sapply(snps.chr1$CHR_POS, get.bin)

get.CC.score <- function(x){
  return(CCscore[x,"score"])
}
snps.chr1$CC.SCORE <- sapply(snps.chr1$BIN, get.CC.score)

get.hist.score <- function(x){
  return(bins.100kb[x,"score"])
}
snps.chr1$HIST.SCORE <- sapply(snps.chr1$BIN, get.hist.score)

plot.CC <- plot.CC + geom_point(data=snps.chr1, aes(x=BIN,y=CC.SCORE))
plot.hist <- plot.hist + geom_point(data=snps.chr1, aes(x=BIN,y=HIST.SCORE))

title <- ggdraw() + 
  draw_label(
    "Chromosome 1 (CLL_110)",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

plot_grid(title, plot.CC, plot.hist, nrow = 3, hjust=0, align="v", rel_heights = c(0.1,1,1))



#### PROVES

