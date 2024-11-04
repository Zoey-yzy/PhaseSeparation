"map protein to ensemble gene name"
library("AnnotationDbi")
library("org.Hs.eg.db")
library("readxl")
data <- read_excel("/Users/zoey/Desktop/PS/output/idmap.xlsx")
data$ensid = mapIds(org.Hs.eg.db,keys=data$symbol,column="ENSEMBL",keytype="SYMBOL",multiVals="first")
data1 <- data[!(data$ensid==""|is.na(data$ensid)|is.na(data$HGNC)),]
ensid <- data1$ensid
df <- data.frame(ensid)
write.csv(df,"/Users/zoey/Desktop/PS/output/ensid.csv")
