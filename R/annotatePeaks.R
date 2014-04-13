annotatePeaks <- function(inputfile, genome, up=NULL, down=NULL){
      ## Check up and down
      if(is.null(up) | is.null(down)) stop("Missing up or down args.")
      
      path <- system.file("extdata",package="GSCA")
      load(paste0(path,"/allreffile.rda"))
      
      ### Read in genome file
      if(genome=="mm8") {
            reffile <- allreffile$mm8
      } else if(genome=="hg18"){
            reffile <- allreffile$hg18
      } else if(genome=="mm9"){
            reffile <- allreffile$mm9
      } else if(genome=="hg19"){
            reffile <- allreffile$hg19
      } else {
            stop("Cannot find genome. Please give 'mm8','mm9','hg18' or 'hg19'")
      }
      
      if(nrow(inputfile)>0){
            
            EntrezGeneID <- rep("0",nrow(inputfile))
            
            ### Change chr notation to match
            chr <- inputfile[,1]
            if(genome=="hg18" | genome=="hg19"){
                  chr[chr=="chrX"] <- "chr23"
                  chr[chr=="chrY"] <- "chr24"
            } else if(genome=="mm8" | genome=="mm9"){
                  chr[chr=="chrX"] <- "chr20"
                  chr[chr=="chrY"] <- "chr21"
            }
            
            ### Adjust reffile according to up-down TSS limits
            TSSlow <- reffile[,6] - up
            TSShigh <- reffile[,6] + down
            TSSlow[reffile[,5]=="-"] <- reffile[reffile[,5]=="-",7] - down
            TSShigh[reffile[,5]=="-"] <- reffile[reffile[,5]=="-",7] + up
            ann <- cbind(reffile[,1:4],TSSlow,TSShigh)
            colnames(ann) <- c("EntrezGeneID","GeneName","NCBIRefseqID","Chr",
                               "TSSlow","TSShigh")
            
            ### Find overlaps
            for(i in 1:nrow(inputfile)){
                  tmp <- ann[strsplit(chr[i],"chr")[[1]][2]==ann[,"Chr"],]
                  
                  start <- as.numeric(inputfile[i,2])
                  end <- as.numeric(inputfile[i,3])
                  index <- ((end <= tmp[,"TSShigh"] & end >= tmp[,"TSSlow"]) |
                                  (start <= tmp[,"TSShigh"] & start >= tmp[,"TSSlow"]) |
                                  (start <= tmp[,"TSSlow"] & end >= tmp[,"TSShigh"]))
                  if(sum(index) > 0){
                        EntrezGeneID[i] <- paste(unique(tmp[index,1]),collapse=",")
                  } else {
                        EntrezGeneID[i] <- "-9"
                  }
            }}
      
      inputfile <- cbind(inputfile,EntrezGeneID)
}
