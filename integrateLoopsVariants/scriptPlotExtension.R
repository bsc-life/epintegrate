# LIBRARIES

library(ggplot2)
library(dplyr)


# LOAD DATA

global <- read.delim("~/Documents/Scripts/Analysis/New/Results/LoopsVariants.csv")


# PLOTS

# Number of loops vs. kb extended in loops regions

global.t <- as.data.frame(table(global$KB.EXTENDED))
ggplot(global.t, aes(x=Var1, y=Freq, group = 1)) + geom_point() + geom_line() + labs(x="Kb extended")

# Number of unique SNPs vs. kb extended in loops regions

x <- c(0,10,20,30,40,50)
y <- c()

for (ext in x){
  s <- filter(global, KB.EXTENDED == as.integer(ext))
  us <- unique(s[,7:8])
  y <- c(y,dim(us)[1])
}

df <- data.frame(x,y)

ggplot(df, aes(x,y, group = 1)) + geom_point() + geom_line() + labs(x="Kb extended", y="SNPS")

# Number of unique SNPs from GWAS catalog database vs. kb extended in loops regions

y1 <- c()

for (ext in x){
  s <- filter(global, KB.EXTENDED == as.integer(ext))
  gwas <- dim(unique(filter(s, DB == "GWAS")[,7:8]))[1]
  y1 <- c(y1,gwas/dim(unique(s[,7:8]))[1])
}

df1 <- data.frame(x,y1)
ggplot(df1, aes(x,y1, group = 1)) + geom_point() + geom_line() + labs(x="Kb extended", y="SNPS GWAS / SNPS Total")

# Number of unique SNPs from GWAS catalog databse with significant p-value ( <= 8*10**(-8)) vs. kb extended in loops regions

y2 <- c()

for (ext in x){
  s2 <- filter(global, KB.EXTENDED == as.integer(ext))
  sig <- filter(s2,  DB == "GWAS" & P.VALUE <= 8*(10**(-8)))
  y2 <- c(y2,dim(unique(sig[,7:8]))[1]/dim(unique(s2[,7:8]))[1])
}

df2 <- data.frame(x,y2)
ggplot(df2, aes(x,y2, group = 1)) + geom_point() + geom_line() + labs(x="Kb extended", y="SNPS significant / SNPS Total")
