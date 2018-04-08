library(GEOquery)
citation("GEOquery")
gds <- getGEO("GDS5260", GSEMatrix = TRUE)
#Description:
head(Meta(gds))
Table(gds)[1:5,]
Columns(gds)[,1:3]
#convert to ExporessionSet (might be easier to work with)
eset <- GDS2eSet(gds,do.log2=TRUE)
#convert to MAList
gpl <- getGEO(Meta(gds)$platform)
MA <- GDS2MA(gds,GPL=gpl)
class(MA)

# GSE
gse <- getGEO("GSE47598", GSEMatrix = FALSE)
head(Meta(gse))
# names of all the GSM objects contained in the GSE
names(GSMList(gse))
# and get the first GSM object on the list
GSMList(gse)[[1]]
# and the names of the GPLs represented
names(GPLList(gse))
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL10558'},GSMList(gse))
length(gsmlist)
Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]
probesets <- Table(GPLList(gse)[[1]])$ID
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
      {tab <- Table(x)
      mymatch <- match(probesets,tab$ID_REF)
      return(tab$VALUE[mymatch])
      }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]
