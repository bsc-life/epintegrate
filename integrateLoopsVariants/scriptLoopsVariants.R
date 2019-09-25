setwd("~/Documents/Scripts/Analysis/New")

# LIBRARIES

library("dplyr")
library("SparseM")
library("plsgenomics")
library("gplots")
library("pheatmap") 
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library(viridis)



# PATHS

# Loops from Homer and Hiccups
path <- getwd()
path.loops.homer <- "./Loops_Homer"
path.loops.hiccups <- "./Loops_Hiccups"

# Create directories to store the results
dir.create("Results")
dir.create("./Results/ResultsHomer")
dir.create("./Results/ResultsHiccups")

# Paths to results directories
path.results.hiccups <- "./Results/ResultsHiccups"
path.results.homer <- "./Results/ResultsHomer"



# LOAD GWAS

#gwas.data <- read.delim("~/Documents/Analysis_GWAS/gwas-association-downloaded_2019-07-04-EFO_0000095-withChildTraits.tsv", na = c("","NA","NR"), stringsAsFactors = FALSE)
#gwas.data <- as.data.frame(gwas.data)

# Filtering repeated reported SNPS, taking the ones with lowest p-values
#new <- data.frame(stringsAsFactors = FALSE)
#for (i in unique(gwas.data$SNPS)){
#  ss <- filter(gwas.data, SNPS == i)
#  m <- min(ss$P.VALUE)
#  new <- rbind(new, filter(ss, P.VALUE == m))
#}

#gwas.data <- new
#gwas.data <- gwas.data[with(gwas.data, order(gwas.data$CHR_ID, gwas.data$CHR_POS)),]

#write.table(gwas.data, "gwas.data.csv", sep="\t", row.names = FALSE)


### NEW !!! ###
gwas.data <- read.delim("~/Documents/Scripts/Analysis/New/variants.data.csv", na = c("","NA","NR"), stringsAsFactors = FALSE)


# FIND INTERSECTION LOOPS AND VARIANTS

files.homer <- list.files(path.loops.homer,full.names = TRUE)
names.homer <- list.files(path.loops.homer)

files.hiccups <- list.files(path.loops.hiccups,full.names = TRUE)
names.hiccups <- list.files(path.loops.hiccups)

# Function to intersect loops with variants, specifying cell line, tool, resolution of the loop and KB extended
find.intersect <- function(fin, gwas.data, CELL, TOOL, RESOLUTION, KB.EXTENDED){
  hic.loops <- read.delim(fin, header=FALSE, comment="#", stringsAsFactors = FALSE)
  hic.loops <- as.data.frame(hic.loops)
  
  intxs <- data.frame(stringsAsFactors = FALSE)
  
  for(i in 1:nrow(gwas.data)) {
    vc <- gwas.data[i,"CHR_ID"]
    vp <- gwas.data[i,"CHR_POS"]
    ext <- KB.EXTENDED*1000
    s <- filter(hic.loops, ((vc == V1 & V2-ext <= vp & V3+ext >= vp)) | ((vc == V4 & V5-ext <= vp & V6+ext >= vp))) 
    if (dim(s)[1] > 0)
      intxs <- unique(rbind(intxs,data.frame(s[1:6],gwas.data[i,],CELL,TOOL,RESOLUTION,KB.EXTENDED,row.names = NULL)))
  }
  
  return(intxs)
}

# KBs to extend at each end of the chromatin contacts (loops) regions
exts <- c(0,10,20,30,40,50)

# Homer loops
outdir <- "./Results/ResultsHomer"

# Run the function for all Homer loops and write the loops with a variant in a new table
for (i in 1:length(files.homer)){
  
  fin <- files.homer[i]
  name.complete <- unlist(strsplit(names.homer[i],"[.]"))
  nm <- name.complete[1] #obtain cell line from file name
  res <- name.complete[2] #obtain resolution of the loop from file name
  
  for (j in 1:length(exts)){ #for the different extensions of the loops regions
    intxs <- find.intersect(fin,gwas.data,nm,"homer",res,exts[j]) #find loops with variants
    if (dim(intxs)[1] > 0){
      if (file.exists(paste(outdir, paste(nm,"csv",sep = "."), sep = "/"))){
        write.table(intxs,paste(outdir, paste(nm,"csv",sep = "."), sep = "/"), sep="\t", col.names = FALSE, row.names = FALSE, append = TRUE)
      } else {
        write.table(intxs,paste(outdir, paste(nm,"csv",sep = "."), sep = "/"), sep="\t", col.names = TRUE, row.names = FALSE, append = TRUE)
      }
    }
  }
  
}

# HiCCUPS loops
outdir <- "./Results/ResultsHiccups"

# Run the function for all HiCCUPS loops and write the loops with a variant in a new table
for (i in 1:length(files.hiccups)){
  
  fin <- files.hiccups[i]
  name.complete <- unlist(strsplit(names.hiccups[i],"[.]"))
  nm <- name.complete[1] #obtain cell line from file name
  res <- "-"
  
  for (j in 1:length(exts)){ #for the different extensions of the loops regions
    intxs <- find.intersect(fin,gwas.data,nm,"hiccups",res,exts[j]) #find loops with variants
    if (count(intxs) > 0){
      if (file.exists(paste(outdir, paste(nm,"csv",sep = "."), sep = "/"))){
        write.table(intxs,paste(outdir, paste(nm,"csv",sep = "."), sep = "/"), sep="\t", col.names = FALSE, row.names = FALSE, append = TRUE)
      } else {
        write.table(intxs,paste(outdir, paste(nm,"csv",sep = "."), sep = "/"), sep="\t", col.names = TRUE, row.names = FALSE, append = TRUE)
      }
    }
  }
  
}



# GLOBAL DATA FRAME
# Join all the results (intersecting loops with variants) into a single file with all the information

files.results.hiccups <- list.files(path.results.hiccups)
files.results.homer <- list.files(path.results.homer)

global <- data.frame(stringsAsFactors = FALSE)

append.data <- function(path,file,db){
  complete.path <- paste(path,file,sep="/")
  dat <- read.delim(complete.path)
  db <- rbind(db,dat)
  return(db)
}

for (f1 in files.results.hiccups){
  global <- append.data(path.results.hiccups,f1,global)
}

for (f2 in files.results.homer){
  global <- append.data(path.results.homer,f2,global)
}

global <- unique(global)

write.table(global, "Results/LoopsVariants.csv", sep="\t", row.names = FALSE)



