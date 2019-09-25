# LIBRARIES
library(dplyr)

# PATHS
path <- getwd()
path.loops.homer <- "./Loops_Homer"
path.loops.hiccups <- "./Loops_Hiccups"

# LOAD DATA
files.homer <- list.files(path.loops.homer,full.names = TRUE)
names.homer <- list.files(path.loops.homer)
files.hiccups <- list.files(path.loops.hiccups,full.names = TRUE)
names.hiccups <- list.files(path.loops.hiccups)
files.all <- c(files.homer,files.hiccups)
names.all <- c(names.homer, names.hiccups)

all.loops <- data.frame(stringsAsFactors = FALSE)

for (i in 1:length(files.all)){
  loops <- read.delim(files.all[1], header=FALSE, comment="#", stringsAsFactors = FALSE)
  cell <- as.character(unlist(strsplit(names.all[i],"[.]"))[1])
  all.loops <- rbind(all.loops, data.frame(dplyr::select(loops, num_range("V", 1:6)),cell))
}

# WRITE FILE
write.table(all.loops, "./AllLoops.csv", sep="\t", row.names = FALSE)
