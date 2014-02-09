tabSearch <- function(keyword,chipdata,option="OR"){
    if(chipdata == "hgu133a") {
        if (!require(Affyhgu133aExpr)) {
            stop("Affyhgu133aExpr Package is not found")
        } else {
            data(Affyhgu133aExprtab)
            tab <- Affyhgu133aExprtab
        }
    } else if(chipdata == "moe430"){
          if (!require(Affymoe430Expr)) {
                stop("Affymoe430Expr Package is not found")
          } else {
                data(Affymoe430Exprtab)
                tab <- Affymoe430Exprtab
          }
    } else {
        stop("Please enter valid name for chipdata. Current Supported chipdata: 'hgu133a' and 'moe430'")
    }

    if(is.null(keyword)) stop("Please enter keyword.")
    keyword <- as.character(keyword)
    keyword <- unique(c(keyword,unlist(strsplit(keyword,";"))))

    if(option=="OR"){
        tmp <- NULL
        for(i in 1:length(keyword)){
            tmp <- rbind(tmp,tab[unique(c(grep(keyword[i],tab$SampleType,ignore.case=T),
                                 grep(keyword[i],tab$ExperimentID,ignore.case=T))),])
        }
    } else if(option=="AND"){
        tmp <- NULL
        for(i in 1:length(keyword)){
            tmp <- c(tmp,unique(c(grep(keyword[i],tab$SampleType,ignore.case=T),
                           grep(keyword[i],tab$ExperimentID,ignore.case=T))))

        }
        tmp <- table(tmp)
        tmp <- tab[as.numeric(names(tmp)[tmp == length(keyword)]),]
    } else {
        stop("Please enter 'OR' or 'AND' for option.")
    }

    if(nrow(tmp)==0) stop("No matching samples found for keyword.")
    tmpdex <- unique(tmp$SampleType)

    output1 <- rep("0",length(tmpdex))
    output2 <- output1
    output3 <- rep(0,length(tmpdex))
    for(i in 1:length(tmpdex)){
        output2[i] <- tmpdex[i]
        output1[i] <- paste(unique(tmp[which(tmp$SampleType==tmpdex[i]),3]),
                             collapse=";")
        output3[i] <- sum(tab$SampleType==tmpdex[i])
    }
    output <- data.frame(output1,output2,output3,
                         stringsAsFactors=FALSE)
    colnames(output) <- c("ExperimentID","SampleType","SampleCount")
    return(output)
}
