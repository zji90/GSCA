ConstructTG <- function(annonPeaksOut,limmaOut){
      
      chip <- unique(unlist(strsplit(as.character(
            annonPeaksOut[,"EntrezGeneID"]),",")))
      chip <- chip[!(chip %in% c("-1","-9"))]
      
      ge <- limmaOut[!is.na(limmaOut[,1]),]
      upge <- ge[ge[,"t"] > 0,1]
      dnge <- ge[ge[,"t"] < 0,1]
      upge2 <- unique(upge[!(upge %in% dnge)])
      dnge2 <- unique(dnge[!(dnge %in% upge)])
      
      uptg <- as.character(upge2[upge2 %in% chip])
      dntg <- as.character(dnge2[dnge2 %in% chip])
      
      return(list(PosTG=uptg,NegTG=dntg))
}
