# LIBRARIES
library(dplyr)

# PATHS
path <- getwd()
path.loops.homer <- "./Loops_Homer"
path.loops.hiccups <- "./Loops_Hiccups"

# LOAD DATA
files.homer <- list.files(path.loops.homer,full.names = TRUE)
files.hiccups <- list.files(path.loops.hiccups,full.names = TRUE)
files.all <- c(files.homer,files.hiccups)
    
variants.data <- read.delim("./variants.data.csv", stringsAsFactors = FALSE) # Only SNPs mapped to CLL trait (from GWAS, DisGeNET and ClinVar)
#all.snps <- read.delim("./all_gwas_catalog_extra.tsv", stringsAsFactors = FALSE) #This one has repeated SNPs
all.snps <- read.delim("./all_gwas.data.csv", stringsAsFactors = FALSE) #For repeated SNPs: keeping the one with lowest p-value

all.loops <- data.frame(stringsAsFactors = FALSE)

for (fin in files.all){
  loops <- read.delim(fin, header=FALSE, comment="#", stringsAsFactors = FALSE)
  all.loops <- rbind(all.loops, dplyr::select(loops, num_range("V", 1:6)))
}

all.snps.unique <- unique(dplyr::select(all.snps, CHR_ID, CHR_POS)) ### Unique or all???

set.seed(5)
stats <- c()

for (i in 1:3){
  x <- sample(1:nrow(all.snps.unique), 244, replace = FALSE)
  intxs <- 0
  for(i in 1:length(x)) {
    vc <- all.snps.unique[x[i],"CHR_ID"]
    vp <- as.numeric(all.snps.unique[x[i],"CHR_POS"])
    s <- filter(all.loops, ((vc == V1 & V2 <= vp & V3 >= vp)) | ((vc == V4 & V5 <= vp & V6 >= vp)))
    intxs <- intxs + dim(s)[1]
  }
  stats <- c(stats,intxs)
}

stats <- as.data.frame(stats)
write.table(stats, "statsLoopsAllGWAS.csv", sep = "\t", quote = FALSE, row.names = FALSE)
