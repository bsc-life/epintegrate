# LIBRARIES

library(dplyr)
library(plyr)


# LOAD DATA

# From GWAS catalog
gwas.data <- read.delim("~/Documents/Scripts/Analysis/New/gwas.data.csv", stringsAsFactors = FALSE) # Already filtered

# From DisGeNET
dgn.data <- read.delim("~/Documents/DataVariants/C0023434_disease_vda_summary.tsv", stringsAsFactors = FALSE)

# From ClinVar
clinvar.data <- read.delim("~/Documents/DataVariants/clinvar_result.txt", stringsAsFactors=FALSE)
clinvar.data <- filter(clinvar.data, grepl('Chronic lymphocytic leukemia', Condition.s.))


# SELECT COLUMNS

gwas.data.s <- select(gwas.data, CHR_ID, CHR_POS, REPORTED.GENE.S., MAPPED_GENE, SNPS, CONTEXT, RISK.ALLELE.FREQUENCY, P.VALUE, OR.or.BETA)
dgn.data.s <- select(dgn.data, Variant, Gene, Chr, Position, Consequence)
clinvar.data.s <- select(clinvar.data, Gene.s., GRCh38Chromosome, GRCh38Location, VariationID)

colnames(dgn.data.s) <- c("SNPS", "MAPPED_GENE", "CHR_ID", "CHR_POS", "CONTEXT")
colnames(clinvar.data.s) <- c("MAPPED_GENE", "CHR_ID", "CHR_POS", "SNPS")

gwas.data.s$CHR_ID <- as.character(gwas.data.s$CHR_ID)
dgn.data.s$CHR_ID <- as.character(dgn.data.s$CHR_ID)
clinvar.data.s$CHR_ID <- as.character(clinvar.data.s$CHR_ID)

gwas.data.s$CHR_POS <- as.numeric(gwas.data.s$CHR_POS)
dgn.data.s$CHR_POS <- as.numeric(dgn.data.s$CHR_POS)
clinvar.data.s$CHR_POS <- as.numeric(clinvar.data.s$CHR_POS)

gwas.data.s$DB <- "GWAS"
dgn.data.s$DB <- "DisGeNET"
clinvar.data.s$DB <- "ClinVar"


# ALL VARIANTS DATAFRAME

all.variants <- data.frame(stringsAsFactors = FALSE)
all.variants <- rbind(all.variants, gwas.data.s)

# Checking
# ints <- intersect(all.variants[1:2],dgn.data.s[3:4]) 
# 106
# 228 - 106 = 122
# 119 + 122 = 241

for (i in 1:nrow(dgn.data.s)){
  chr_pos <- dgn.data.s[i,"CHR_POS"]
  chr_id <- dgn.data.s[i,"CHR_ID"]
  s <- filter(all.variants, CHR_ID == chr_id & CHR_POS == chr_pos)
  if (dim(s)[1] == 0) # variant not already in all.variants
    all.variants <- rbind.fill(all.variants,dgn.data.s[i,])
}

for (j in 1:nrow(clinvar.data.s)){
  chr_pos <- clinvar.data.s[j,"CHR_POS"]
  chr_id <- clinvar.data.s[j,"CHR_ID"]
  s <- filter(all.variants, CHR_ID == chr_id & CHR_POS == chr_pos)
  if (dim(s)[1] == 0) # variant not already in all.variants
    all.variants <- rbind.fill(all.variants,clinvar.data.s[j,])
}


# CONTEXT RENAME TO UNIFORMIZE

rename <- function(x){
  x <- paste(unlist(strsplit(x, " +")), collapse="_")
  return(x)
}
all.variants$CONTEXT <- sapply(all.variants$CONTEXT, rename)


# REORDER BY CHROMOSOMIC LOCATION

f1 <- filter(all.variants, CHR_ID != "X" & CHR_ID != "Y" & CHR_ID != "MT")
f1$CHR_ID <- as.numeric(f1$CHR_ID)
f1 <- f1[with(f1, order(f1$CHR_ID,f1$CHR_POS)),]
f1$CHR_ID <- as.character(f1$CHR_ID)

f2 <- filter(all.variants, CHR_ID == "X")
f2 <- f2[with(f2, order(f2$CHR_POS)),]

f3 <- filter(all.variants, CHR_ID == "Y")
f3 <- f3[with(f3, order(f3$CHR_POS)),]

f4 <- filter(all.variants, CHR_ID == "MT")
f4 <- f4[with(f4, order(f4$CHR_POS)),]

all.variants <- rbind(f1,f2,f3,f4)


# SAVE ALL VARIANTS IN A UNIQUE FILE

write.table(all.variants, "variants.data.csv", sep="\t", row.names = FALSE)


