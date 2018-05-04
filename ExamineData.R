# Install GEOquery library:
# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")
library(GEOquery)
library("preprocessCore")
citation("GEOquery")
# GSE
gse <- getGEO("GSE47598", GSEMatrix = FALSE)
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL10558'},GSMList(gse))
PvalTable <- Table(gsmlist[[1]])
Columns(gsmlist[[1]])[1:5,]
probesets <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
      {tab <- Table(x)
      mymatch <- match(probesets,tab$ID_REF)
      return(tab$VALUE[mymatch])
      }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix[is.na(data.matrix)] = 0

#Data Processing mentioned in the Paper
#data.matrix = data.matrix[rowSums(data.matrix == 0) <= 3, ]
data.matrix <- log2(data.matrix)
data.matrix[is.na(data.matrix)] = 0
data.matrix = normalize.quantiles(data.matrix,copy=TRUE)

  
data.labels <- lapply(GSMList(gse),function(x) {Meta(x)$characteristics_ch1}) == "group: Psoriasis Patient"
data.rowNames <- probesets 
data.Pval <- rep(1, times = nrow(data.matrix))
for (ii in range(0, nrow(data.matrix))) {
  data.Pval[ii] = PvalTable[PvalTable[, 1]== data.rowNames[ii], 3]
}

# Recreating results
#data.Filmatrix = data.matrix[data.Pval <= 0.05, ]
#data.Filmatrix = data.Filmatrix[rowSums(data.Filmatrix == 0) <= 3, ]
infected <- data.matrix[, data.labels]
ninfected <- data.matrix[, !data.labels]
meanInfected <- rowMeans(infected)
meanNinfected <- rowMeans(ninfected)
expressChange <- (meanInfected )/meanNinfected
#List of Upregulated genes 
upregSort <- order(expressChange, decreasing = TRUE)
data.rowNames[upregSort[0:13]] #upregSort maps indexes to sorted indexes
expressChange[upregSort[0:13]]

#testing
leader <- which(data.rowNames =="ILMN_1683678")[[1]]
expressChange[upregSort[leader]]

