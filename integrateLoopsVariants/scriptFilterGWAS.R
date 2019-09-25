input.path <- "~/Documents/Scripts/Analysis/New/all_gwas_catalog_extra.tsv"
output.name <- "all_gwas.data.csv"

# LIBRARIES
library(dplyr)

gwas.data <- read.delim(input.path, na = c("","NA","NR"), stringsAsFactors = FALSE)
gwas.data <- as.data.frame(gwas.data)

# Filtering repeated reported SNPS, taking the ones with lowest p-values
new <- data.frame(stringsAsFactors = FALSE)
for (i in unique(gwas.data$SNPS)){
  ss <- filter(gwas.data, SNPS == i)
  m <- min(ss$P.VALUE)
  new <- rbind(new, filter(ss, P.VALUE == m))
}

gwas.data <- new
gwas.data <- gwas.data[with(gwas.data, order(gwas.data$CHR_ID, gwas.data$CHR_POS)),]

write.table(gwas.data, output.name, sep="\t", row.names = FALSE)