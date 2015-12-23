# Analysing ranking of HVG's
source("./functions.R")
source("./DM.R")
library("plyr")

out.dir <- "./HVGresults"
dir.create(out.dir)

for (x in c("Zeisel","Klein")) {
    cur.data <- readRDS(sprintf("%sData.rds", x))
    counts <- cur.data$Counts
    size.facs <- cur.data$SF
    ranking <- data.frame(
                          "DESeq"=character(length=nrow(counts)),
                          "TMM"=character(length=nrow(counts)),
                          "LibSize"=character(length=nrow(counts)),
                          "Deconvolution"=character(length=nrow(counts)),
                          stringsAsFactors=FALSE
                             )
    for (y in colnames(size.facs)) {
             norm.counts <- as.data.frame(t(t(counts) / size.facs[,y]))
             features <- featureCalc(norm.counts, htseq=FALSE, 1)
             dm <- DM(meanGenes=features$mean,CV2Genes=features$cv2)
             ranking[,y] <- rownames(features[order(-dm),]) 
             }

    for (i in c(500,1000,2000)) {
        comparisonMatrix <- compareHVG(ranking,i)
        write.table(comparisonMatrix,file.path(out.dir,paste(x,"HVGTop",i,".txt",sep="")))
        }

}


