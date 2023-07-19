ulm.ebayes <- function(signal, matchMatrix, use.intercept=TRUE, minSize=5){
    signal <- signal[names(signal) %in% row.names(matchMatrix)]
    matchMatrix <- matchMatrix[names(signal),]
    matchMatrix <- matchMatrix[,which(colSums(matchMatrix!=0L)>=minSize)]
    if(use.intercept){
      modelmtx <- model.matrix(~signal)
    }else{
      modelmtx <- model.matrix(~0+signal)
    }
    lmfit <- lmFit(as.matrix(t(matchMatrix)), modelmtx)
    ulmres <- as.data.frame(topTable(eBayes(lmfit), coef="signal", number = Inf)[,c(1,4,5)])
    colnames(ulmres) <- c("score", "p", "padj")
    ulmres
  }
 
runulm <- function(DARmat,
                   matchMtx){  
   
  ptm <- proc.time()
  ulm <- ulm.ebayes(DARmat, matchMtx)
  runtime <- proc.time()-ptm
  
  return(list(ulm, runtime))
}