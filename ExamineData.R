library(GEOquery)
citation("GEOquery")
# GSE
gse <- getGEO("GSE47598", GSEMatrix = FALSE)
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL10558'},GSMList(gse))
Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]
probesets <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
      {tab <- Table(x)
      mymatch <- match(probesets,tab$ID_REF)
      return(tab$VALUE[mymatch])
      }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
#data.matrix <- log2(data.matrix)
data.matrix[is.na(data.matrix)] = 0
data.matrix[5:12,]
data.labels <- lapply(GSMList(gse),function(x) {Meta(x)$characteristics_ch1}) == "group: Psoriasis Patient"
data.rowNames <- probesets

# Recreating results
infected <- data.matrix[, data.labels]
ninfected <- data.matrix[, !data.labels]
meanInfected <- rowMeans(infected)
meanNinfected <- rowMeans(ninfected)
expressChange <- (meanInfected - meanNinfected)/meanNinfected
#List of Upregulated genes
upregSort <- order(expressChange, decreasing = TRUE)
data.rowNames[upregSort[0:13]]
expressChange[upregSort[0:13]]

#testing
which(data.rowNames =="ILMN_1683678")[[1]]

