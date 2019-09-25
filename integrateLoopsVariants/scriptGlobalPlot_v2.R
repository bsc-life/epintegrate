# LIBRARIES

library("dplyr")
library("tidyr")
library("ggplot2")
library("viridis")
library("gridExtra")
library("cowplot")

# LOAD DATA

global <- read.delim("~/Documents/Scripts/Analysis/New/Results/LoopsVariants.csv", stringsAsFactors = FALSE)
gwas.data <- read.delim("~/Documents/Scripts/Analysis/New/variants.data.csv", stringsAsFactors = FALSE) #244 variants from different databases
complementary.data <- read.delim("~/Documents/Scripts/Analysis/New/annotationComplementaryGCVariants.csv", stringsAsFactors = FALSE)

gwas.data <- dplyr::select(gwas.data, CHR_ID, CHR_POS, REPORTED.GENE.S., MAPPED_GENE, SNPS, CONTEXT, RISK.ALLELE.FREQUENCY, P.VALUE, OR.or.BETA, DB)

# MATRIX

get.matrix <- function(gwas.data, total.data){
  cell.lines <- sort(unique(total.data$CELL))
  mat <- matrix(0,nrow=length(cell.lines),ncol=nrow(gwas.data))
  for (i in 1:length(cell.lines)){
    for (j in 1:nrow(gwas.data)){
      snps <- gwas.data[j,"SNPS"]
      chr_id <- gwas.data[j,"CHR_ID"]
      chr_pos <- gwas.data[j,"CHR_POS"]
      s <- unique(filter(total.data, CELL==cell.lines[i], SNPS==snps, CHR_ID==chr_id, CHR_POS==chr_pos))
      if (dim(s)[1] > 0)
        mat[i,j] = as.numeric(dim(s)[1])
    }
  }
  return(mat)
}

global.0 <- filter(global, KB.EXTENDED == 0)
mat.0 <- get.matrix(gwas.data,global.0)

get.matrix.complementary <- function(variants.data, complementary.data){
  cell.lines <- sort(unique(complementary.data$CELL))
  mat <- matrix(0,nrow=length(cell.lines),ncol=nrow(variants.data))
  for (i in 1:length(cell.lines)){
    for (j in 1:nrow(variants.data)){
      snps <- variants.data[j,"SNPS"]
      chr_id <- variants.data[j,"CHR_ID"]
      chr_pos <- variants.data[j,"CHR_POS"]
      s <- unique(filter(complementary.data, CELL==cell.lines[i], SNPS==snps, CHR_ID==chr_id, CHR_POS==chr_pos))
      if (dim(s)[1] > 0)
        mat[i,j] = s[1,"gc"]
    }
  }
  return(mat)
}

complementary.matrix <- get.matrix.complementary(variants.data = gwas.data, complementary.data = complementary.data)

# PLOT

base_size = 9

colours <- c("intergenic_variant"="#F8766D","3_prime_UTR_variant"="#E69F00","regulatory_region_variant"="#00BFC4","intron_variant"="#56B4E9","synonymous_variant"="#00BA38","non_coding_transcript_exon_variant"="#F564E3",
             "5_prime_UTR_variant"="#F0E442","splice_region_variant"="#619CFF","missense_variant"="#CC79A7","TF_binding_site_variant"="#9999CC")

gwas.data$name <- paste("chr",gwas.data$CHR_ID," ",gwas.data$CHR_POS, sep = "")
gwas.data$name <- factor(gwas.data$name, levels=unique(gwas.data$name))

gwas.data$SNPS <- factor(gwas.data$SNPS, levels = unique(gwas.data$SNPS))

currentMat <- mat.0 # !!!

df <- as.data.frame(currentMat)
colnames(df) <- paste("chr",gwas.data$CHR_ID," ",gwas.data$CHR_POS, sep = "")
df$CELL <- sort(unique(global$CELL))
df.g <- gather(df,SNP,Count,-CELL)
df.g$SNP <- factor(df.g$SNP, levels = unique(df.g$SNP))

gg <- ggplot(df.g, aes(x=SNP,y=CELL,fill=Count)) + geom_tile(colour="white") 
gg <- gg + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_fill_gradient(low = "white", high = "steelblue")

df.compl <- as.data.frame(complementary.matrix)
colnames(df.compl) <- paste("chr",gwas.data$CHR_ID," ",gwas.data$CHR_POS, sep = "")
df.compl$CELL <- sort(unique(complementary.data$CELL))
df.compl.g <- gather(df.compl,SNP,ComplementaryGC,-CELL)
df.compl.g$SNP <- factor(df.compl.g$SNP, levels = unique(df.compl.g$SNP))

colours.compl <- c("0"="white","exon"="#F564E3","intron"="#00BFC4","1to5kb"="#E69F00","intergenic"="#00BA38")

cc <- ggplot(df.compl.g, aes(x=SNP,y=CELL,fill=ComplementaryGC)) + geom_tile(colour="white") 
cc <- cc + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + scale_fill_manual(values=colours.compl) 

# Cut off p-value: 8*(8**(-10))
# log2(8*(8**(-10))) = -27
pv <- ggplot(gwas.data, aes(x = name, y = 0.5, fill = log2(P.VALUE))) + geom_tile(colour="white") 
pv <- pv + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(angle = 0, margin=margin(0,0,0,0), vjust = .5), axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_fill_viridis(limits = c(-200, -27), na.value="grey80", direction = -1) + ylab("P-value")

or <- ggplot(gwas.data, aes(x = name, y = 0.5, fill = log2(OR.or.BETA))) + geom_tile(colour="white") 
or <- or + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(angle = 0, margin=margin(0,0,0,0), vjust = .5), axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_fill_viridis(na.value="grey80")  + ylab("OR")

raf <- ggplot(gwas.data, aes(x = name, y = 0.5, fill = RISK.ALLELE.FREQUENCY)) + geom_tile(colour="white") 
raf <- raf + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(angle = 0, margin=margin(0,0,0,0), vjust = .5), axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_fill_viridis(na.value="grey80") + ylab("RAF")

db <- ggplot(gwas.data, aes(x = name, y = 0.5, fill = DB)) + geom_tile(colour="white")
db <- db + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(angle = 0,margin=margin(0,0,0,0), vjust = .5), axis.text.y = element_blank(), axis.ticks = element_blank()) + ylab("DB")

gc <- ggplot(gwas.data, aes(x = name, y = 0.5, fill = CONTEXT)) + geom_tile(colour="white") 
gc <- gc + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 90, hjust = 1, vjust = .5, colour = "grey50"), axis.title.y = element_text(angle = 0,margin=margin(0,0,0,0), vjust = .5), axis.text.y = element_blank(), axis.ticks = element_blank()) + ylab("Genomic \n context") #+ scale_fill_manual(values=colours) 

#Format legends
gg <- gg + labs(fill="Counts") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
cc <- cc + labs(fill="Complementary GC") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
pv <- pv + labs(fill="Log P-value") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
or <- or + labs(fill="Log OR") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
raf <- raf + labs(fill="RAF") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
db <- db + labs(fill="DB") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.7, "cm"))
gc <- gc + labs(fill="Genomic context") + theme(legend.title = element_text(size = 10), legend.text = element_text(size = 7), legend.key.size = unit(.3, "cm"))

#Get legends separately from plots
l.gg <- cowplot::get_legend(gg)
l.cc <- cowplot::get_legend(cc)
l.pv <- cowplot::get_legend(pv)
l.or <- cowplot::get_legend(or)
l.raf <- cowplot::get_legend(raf)
l.db <- cowplot::get_legend(db)
l.gc <- cowplot::get_legend(gc)

#Remove legends from plots
gg <- gg + theme(legend.position = "none")
cc <- cc + theme(legend.position = "none")
pv <- pv + theme(legend.position = "none")
or <- or + theme(legend.position = "none")
raf <- raf + theme(legend.position = "none")
db <- db + theme(legend.position = "none")
gc <- gc + theme(legend.position = "none")


all_plots <- plot_grid(gg, cc, pv, or, raf, db, gc, hjust=0, nrow = 7, rel_heights = c(3/2, 3/2, 2/7, 2/7, 2/7, 2/7, 6/7), align = "v")
all_legends <- plot_grid(l.gg, l.cc, l.pv, l.or, l.raf, l.db, l.gc, nrow = 1, ncol = 7, scale = 1)
plot_grid(all_plots, all_legends, nrow = 2, rel_heights = c(5/7,2/7))


# SAVE PLOT
ggsave('plotGlobal.v2.png', path = '~/Documents/Scripts/Analysis/New/', width = unit(15, "cm"), height = unit(10, "cm"))


