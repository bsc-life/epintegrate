# LIBRARIES
library(dplyr)
library(rlist)

# LOAD DATA
all.loops.data <- read.delim("./AllLoops.csv", stringsAsFactors = FALSE)

# Tots els loops son intra!!
# No hi ha loops al chr MT

# Human Genome Assembly GRCh38
# www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
chr.list <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y") #include MT? No data in fransua/results
dict.length.chr <- c(248956422, 242193529,	198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)

cells <- c("CLL110","CLL12","CLL1525","GCBC","MBC","NBC","PBC")

all.scores <- list()

# SEPARAR PER CELL LINES
for (c in cells){
  scores <- list()
  all.loops <- filter(all.loops.data, cell == "CLL110")
  
  for (chr in 1:length(chr.list)){
    chr.score <- c()
    loops.chr <- filter(all.loops, V1 == chr.list[chr] | V4 == chr.list[chr])
    for (i in seq(1, dict.length.chr[chr], by=100000)){
      
      #Both start and end within the bin
      se.1 <- filter(loops.chr, V2 >= i & V2 < (i + 100000) & V3 >= i & V3 < (i + 100000))
      se.1 <- dplyr::select(se.1, V1, V2, V3)
      colnames(se.1) <- c("CHR","START","END")
      
      se.2 <- filter(loops.chr, V5 >= i & V5 < (i + 100000) & V6 >= i & V6 < (i + 100000))
      se.2 <- dplyr::select(se.2, V4, V5, V6)
      colnames(se.2) <- c("CHR","START","END")
      
      #Only start within the bin
      s.1 <- filter(loops.chr, V2 >= i & V2 < (i + 100000) & V3 >= (i + 100000))
      s.1 <- dplyr::select(s.1, V1, V2, V3)
      if (nrow(s.1) > 0) {
        s.1$V3 <- as.integer(i + 100000)
      }
      colnames(s.1) <- c("CHR","START","END")
      
      s.2 <- filter(loops.chr, V5 >= i & V5 < (i + 100000) & V6 >= (i + 100000))
      s.2 <- dplyr::select(s.2, V4, V5, V6)
      if (nrow(s.2) > 0) {
        s.2$V6 <- as.integer(i + 100000)
      }
      colnames(s.2) <- c("CHR","START","END")
      
      #Only end within the bin
      e.1 <- filter(loops.chr, V2 < i & V3 >= i & V3 < (i + 100000))
      e.1 <- dplyr::select(e.1, V1, V2, V3)
      if (nrow(e.1) > 0) {
        e.1$V2 <- as.integer(i)
      }
      colnames(e.1) <- c("CHR","START","END")
      
      e.2 <- filter(loops.chr, V5 < i & V6 >= i & V6 < (i + 100000))
      e.2 <- dplyr::select(e.2, V4, V5, V6)
      if (nrow(e.2) > 0) {
        e.2$V5 <- as.integer(i)
      } 
      colnames(e.2) <- c("CHR","START","END")
      
      loops.regions <- rbind(se.1,se.2,s.1,s.2,e.1,e.2)
      
      if (i > dict.length.chr[chr]-100000){
        norm <- dict.length.chr[chr] - i
        loops.regions$LEN <- ((loops.regions$END - loops.regions$START)*100000)/norm
      }
      loops.regions$LEN <- loops.regions$END - loops.regions$START
      
      chr.score <- c(chr.score,sum(loops.regions$LEN))
    }
    scores[[chr]] <- chr.score
  }
  all.scores <- list.append(all.scores, scores)
}

# SAVE LIST SCORES
list.save(all.scores, "allScoresCC.RData")

# SAVE SEPARATED BY CHR AND IN .CSV
all.scores <- list.load("allScoresCC.RData")

dir.create("./CCscores")
setwd("./CCscores")

for (i in 1:length(cells)){
  for (j in 1:length(chr.list)){
    l <- all.scores[[i]][[j]]
    df <- data.frame(matrix(unlist(l), byrow=T),stringsAsFactors=FALSE)
    colnames(df) <- c("score")
    write.table(df, paste(cells[i],"chr",chr.list[j],".csv",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
