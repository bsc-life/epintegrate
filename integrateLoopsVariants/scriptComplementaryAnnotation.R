# LOAD LIBRARIES

# BiocManager::install(c("rtracklayer"))
# BiocManager::install(c("annotatr"))
# BiocManager::install("GenomicRanges")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("org.Hs.eg.db")

library(dplyr)
library(rtracklayer)
library(annotatr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(gdata)


# LOAD DATA

LoopsVariants <- read.delim("~/Documents/Scripts/Analysis/New/Results/LoopsVariants.csv")


# GET "COMPLEMENTARY REGIONS" FROM LOOPS - VARIANTS DF

# Annotating the loops regions themselves without extending them at the ends
LoopsVariants.0 <- filter(LoopsVariants, KB.EXTENDED == 0)

# Function to see in which region of the loop is the SNP located
which.region.loop <- function(r){
  if (r["CHR_ID"] == r["V1"] & r["V2"] <= r["CHR_POS"] & r["V3"] >= r["CHR_POS"])
    return("+") # first region of the loop in genomic position order
  else
    return("-") # second region of the loop in genomic position order
}

LoopsVariants.0$LOOP.REGION <- apply(LoopsVariants.0, 1, which.region.loop)

complementary.regions <- data.frame()
snps.regions <- data.frame()

s1 <- filter(LoopsVariants.0, LOOP.REGION == "+")
p1 <- dplyr::select(s1, V1, V2, V3) # SNP located here
c1 <- dplyr::select(s1, V4, V5, V6) # complementary region
colnames(p1) <- c("chrom","chromStart","chromEnd")
colnames(c1) <- c("chrom","chromStart","chromEnd")

s2 <- filter(LoopsVariants.0, LOOP.REGION == "-")
p2 <- dplyr::select(s2, V4, V5, V6) # SNP located here
c2 <- dplyr::select(s2, V1, V2, V3) # complementary region
colnames(p2) <- c("chrom","chromStart","chromEnd")
colnames(c2) <- c("chrom","chromStart","chromEnd")

# Table for loops regions with SNPs
snps.regions <- rbind(snps.regions,unique(p1),unique(p2))

# Table for loops region complementary to SNPs
complementary.regions <- rbind(complementary.regions,unique(c1),unique(c2)) ### unique because KB.EXTENDED

rename <- function(x){
  x <- paste("chr", as.character(x), sep="")
  return(x)
}

snps.regions$chrom <- sapply(snps.regions$chrom, rename)
complementary.regions$chrom <- sapply(complementary.regions$chrom, rename)

# Write table of regions with SNPs
write.table(snps.regions, "snpsRegions.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Write table of complementary regions
write.table(complementary.regions, "complementaryRegions.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Get annotations
an <- build_annotations(genome='hg38', annotations = 'hg38_basicgenes')

# Read complementary regions and annotate
gr.complementary <- read_regions(con = "complementaryRegions.bed", genome = "hg38", format = "bed")
out.complementary <- annotate_regions(regions = gr.complementary, annotations = an)

annotation.complementary.regions <- as.data.frame(out.complementary)
write.table(annotation.complementary.regions, "annotationComplementaryRegions.csv", sep = "\t", quote = FALSE, row.names = FALSE)


# Read loop regions with SNPs and annotate
gr.snps <- read_regions(con = "snpsRegions.bed", genome = "hg38", format = "bed")
out.snps <- annotate_regions(regions = gr.snps, annotations = an)

annotation.snps.regions <- as.data.frame(out.snps)
write.table(annotation.snps.regions, "annotationSNPsRegions.csv", sep = "\t", quote = FALSE, row.names = FALSE)


# Annotate complementary

# Format column with context obtained from annotation package
get.context <- function(x){
  gc <- unlist(strsplit(x, ":", fixed = TRUE))[1]
  return(gc)    
}

annotate.complementary <- function(var,annot){
  
  # Established order
  priorities <- c("exon","5UTR","3UTR","promoter","intron","1to5kb") 
  
  if (var["LOOP.REGION"] == "+") {
    loop.annot <- filter(annot, seqnames == paste("chr", as.character(var["V4"]), sep="") & start == as.integer(var["V5"])+1 & end == as.integer(var["V6"]))
  } else {
    loop.annot <- filter(annot, seqnames == paste("chr", as.character(var["V1"]), sep="") & start == as.integer(var["V2"])+1 & end == as.integer(var["V3"]))
  }
  
  if (dim(loop.annot)[1] > 0){
    loop.annot$gc <- sapply(loop.annot$annot.id, get.context)
    loop.annot$gc <- reorder.factor(loop.annot$gc, new.order = priorities)
    loop.annot <- arrange(loop.annot, gc)
    return(loop.annot[1,])
  } else {
    return(data.frame(gc="intergenic"))
  }
}

loops.vars.compl.annot <- apply(X=LoopsVariants.0, MARGIN=1, FUN=function(var2) annotate.complementary(var2,annot=annotation.complementary.regions))
loops.vars.compl.annot <- rbindlist(loops.vars.compl.annot, fill = TRUE)

annotated.df <- data.frame(LoopsVariants.0,loops.vars.compl.annot)
write.table(annotated.df, "annotationComplementaryGCVariants.csv", sep = "\t", quote = FALSE, row.names = FALSE)

t <- table(annotated.df$gc)
# Exons 0.7, intergenic 0.18, 1to5kb 0.09, introns 0.03


# Check genomic context SNPs

#NOT WORKING WHEN MORE THAN ONE GENE IN MAPPED_GENE OR REPORTED.GENE.S.

check.context <- function(var,annot){
  priorities <- c("exon","5UTR","3UTR","promoter","intron","1to5kb")
  if (var["LOOP.REGION"] == "+") {
    loop.annot <- filter(annot, seqnames == paste("chr", as.character(var["V1"]), sep="") & start == as.integer(var["V2"])+1 & end == as.integer(var["V3"]))
  } else {
    loop.annot <- filter(annot, seqnames == paste("chr", as.character(var["V4"]), sep="") & start == as.integer(var["V5"])+1 & end == as.integer(var["V6"]))
  }
  
  if (dim(loop.annot)[1] > 0){
    s <- filter(loop.annot, as.character(annot.symbol) %in% as.character(var["REPORTED.GENE.S."]) | as.character(annot.symbol) %in% as.character(var["MAPPED_GENE"]))
    #s <- filter(loop.annot, grepl(as.character(annot.symbol),as.character(var["REPORTED.GENE.S."]),fixed=TRUE) | grepl(as.character(annot.symbol),as.character(var["MAPPED_GENE"]),fixed=TRUE))
    gc.check <- (dim(s)[1]>0)
    if (gc.check){
      s$gc <- sapply(s$annot.id, get.context)
      s$gc <- reorder.factor(s$gc, new.order = priorities)
      s <- arrange(s, gc)
      gc.snp <- s[1,"gc"]
      return(data.frame(s[1,],gc.snp=gc.snp,gc.check=gc.check))
    } else {
      loop.annot$gc <- sapply(loop.annot$annot.id, get.context)
      loop.annot$gc <- reorder.factor(loop.annot$gc, new.order = priorities)
      loop.annot <- arrange(loop.annot, gc)
      gc.snp <- loop.annot[1,"gc"]
      return(data.frame(loop.annot[1,],gc.snp=gc.snp,gc.check=gc.check))
    }
  } else {
    gc.snp <- "intergenic"
    gc.check <- (var["CONTEXT"] == "intergenic_variant")
    return(data.frame(gc.snp=gc.snp,gc.check=gc.check))
  }
  
}

loops.vars.gc.check <- apply(X=LoopsVariants.0, MARGIN=1, FUN=function(var2) check.context(var2,annot=annotation.snps.regions))
loops.vars.gc.check <- rbindlist(loops.vars.gc.check, fill = TRUE)

annotated.df.2 <- data.frame(LoopsVariants.0,loops.vars.gc.check)




#####################################

# Enhancers??
# build_annotations(annotations = build_enhancer_annots(genome = "hg38"))
# enh <- build_enhancer_annots(genome = "hg38")

#####################################

# cc <- filter(LoopsVariants, KB.EXTENDED == 0)

# prova3 <- unlist(strsplit(prova2$annot.id, split=":"))

# pp <- LoopsVariants[,7:8]
# pp_count <- unique(pp)

# pppp <- filter(LoopsVariants, KB.EXTENDED == 0)
# pp_count2 <- unique(pppp[,1:8])

# unique.lv <- unique(filter(LoopsVariants, KB.EXTENDED == 0)[,1:8])
# dup.lv <- duplicated(pppp[,1:8])
# unique.lv.counts <- aggregate(list(numdup=rep(1,nrow(pppp))),pppp,length)
# unique.lv.counts

# How many times is duplicated each SNP?
# q <- filter(LoopsVariants, KB.EXTENDED == 0)[,7:8]
# df.counts.snps <- aggregate(list(numdup=rep(1,nrow(q))), q, length)
# df.counts.snps
# max(df.counts.snps$numdup)

#For each cell line
# library(dplyr)
# qq <- filter(LoopsVariants,KB.EXTENDED==0)
# qq.cells <- dplyr::select(qq, CHR_ID, CHR_POS, CELL)
# qq.counts <- unique(qq.cells)

