# Pierre-Luc Germain's TFA method

#' addMotifRelativeWeights
#'
#' Adds a simple relative weight to the motif instances (score divided max other)
#'
#' @param pmoi A GRanges of motif instances, with at least 'score' and 'motif_id' mcol columns.
#' @param minoverlap The minimum number of overlapping nucleotides to consider 
#'   two motif instances as overlapping.
#' @param normalizeScores Logical; whether to normalize scores by the max score for that motif
#'
#' @return The GRanges with an additional 'weight' column
addMotifRelativeWeights <- function(pmoi, minoverlap=3L, normalizeScores=TRUE){
  if(normalizeScores){
    msmax <- max(splitAsList(pmoi$score, pmoi$motif_id))
    pmoi$score <- pmoi$score/msmax[pmoi$motif_id]
  }
  o <- findOverlaps(pmoi, pmoi, minoverlap=minoverlap)
  o <- o[which(from(o)!=to(o) & pmoi$motif_id[from(o)]!=pmoi$motif_id[to(o)])] # remove self-overlaps
  # for each instance that has overlap, identify the max score or other motifs
  pmoi$maxOther <- 1 # when no overlap
  o <- max(splitAsList(pmoi$score[to(o)], from(o))) 
  pmoi$maxOther[as.integer(names(o))] <- pmax(o,1)
  # divide the score of this instance by the maximum other score
  pmoi$weight <- (pmax(pmoi$score,0)/pmoi$maxOther)
  pmoi
}

#' getMotifEnrScores
#' 
#' Get sample-wise motif accessibility scores based on insertion counts flanking
#' competitively-weighted motif matches.
#'
#' @param x A named vector of bam files or fragment files.
#' @param pmoi A GRanges of motif instances, with a 'weight' and a 'motif_id' column.
#' @param wSizes Positive integer of length 2 indicating the window sizes around
#'   motifs. The first number is immediately flanking and counted fully, while
#'   the weight of the second window decays as the inverse of distance
#' @param nf The lower and upper size boundaries of nucleosome-free fragments.
#'   Estimates based on both NF and all fragments will be returned.
#' @param param An optional ScanBamParam object (by default ignores duplicates).
#' @param BPPARAM A BiocParallel param object for multithreading
#' @param scoreFn Function providing the actual scoring
#'
#' @return A summarized experiment if length(x)>1, or a list containing the 
#'   matrices otherwise.
#' @importFrom SummarizedExperiment SummarizedExperiment rowData rowData<-
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom BiocParallel SerialParam bplapply
#' @import S4Vectors IRanges GenomicRanges
#' @export
getMotifEnrScores <- function(x, pmoi, wSizes=c(10,40), nf=c(50,120), param=NULL,
                              BPPARAM=BiocParallel::SerialParam(), scoreFn=.getScores){
  stopifnot(!is.null(pmoi$weight))
  if(is.null(param)){
    reg <- GenomicRanges::reduce(resize(pmoi, width=500, fix="center"))
    param <- ScanBamParam(which=reg,
                          flag=scanBamFlag(isDuplicate=FALSE, isProperPair=TRUE,
                                           isNotPassingQualityControls=FALSE))
  }
  if(length(x)>1){
    # for multiple samples, run separately and aggregate in a SE
    res <- lapply(x, pmoi=pmoi, wSizes=wSizes, nf=nf, param=param,
                  BPPARAM=BPPARAM, FUN=getMotifEnrScores)
    fu <- do.call(cbind, lapply(res, FUN=function(x) x$es[,1]))
    se <- SummarizedExperiment(list(
      NF=do.call(cbind, lapply(res, FUN=function(x) x$nf.es[,1])), # nucleosome-free
      full=do.call(cbind, lapply(res, FUN=function(x) x$es[,1])), # all fragments
      nf.sd=do.call(cbind, lapply(res, FUN=function(x) x$nf.es[,2])), # sd of NF
      full.sd=do.call(cbind, lapply(res, FUN=function(x) x$es[,2])))) # sd of all
    se$nfrags <- sapply(res, FUN=function(x) x$frags) # total frags
    se$nff <- sapply(res, FUN=function(x) x$nff) # total NF frags
    tt <- table(pmoi$motif_id)
    rowData(se)$instances <- as.numeric(tt[row.names(se)])
    return(se)
  }
  message(x)
  if(grepl("\\.bam$",x,ignore.case=TRUE)){
    # read the bam file
    x <- as(readGAlignmentPairs(x, param=param), "GRanges")
  }else{
    # assuming this is a fragment file
    x <- data.table::fread(x, select=1:3)
    x <- GRanges(x$V1, IRanges(x$V2, x$V3))
  }
  ls <- length(x)
  # split by frag size
  message("Read ", ls, " fragments")
  x2 <- x[width(x)>=nf[1] & width(x)<=nf[2]]
  ls2 <- length(x2)
  # get both coverages
  x <- coverage(epiwraps:::.align2cuts(GRanges(x)))[seqlevels(pmoi)]
  x2 <- coverage(epiwraps:::.align2cuts(GRanges(x2)))[seqlevels(pmoi)]
  # run the scoring
  list(
    frags=ls, nff=ls2,
    es=scoreFn(x, pmoi, wSizes=wSizes, BPPARAM=BPPARAM),
    nf.es=scoreFn(x2, pmoi, wSizes=wSizes, BPPARAM=BPPARAM)
  )
}

# Provides weighted insertion counts in windows surrounding motifs
# x is a coverage RleList
# pmoi are the motif instances, including the "weight" mcol column
# wSizes are the size of windows to use (length 2), and wWeight are the weights given to the two sizes 
.getScores <-  function(x, pmoi, wSizes=c(10L,40L), wWeight=c(10,1), minoverlap=3L,
                        BPPARAM=BiocParallel::SerialParam()){
  res <- bplapply(split(pmoi, pmoi$motif_id), BPPARAM=BPPARAM, FUN=function(p){
    message(p$motif_id[1])
    p <- sort(p)
    # count the insertions in the flanking windows (2 sizes) on each side of motif instances
    p$enr.score1 <- (unlist(viewSums(Views(x, flank(p, width=wSizes[1], start=TRUE))))+
                       unlist(viewSums(Views(x, flank(p, width=wSizes[1], start=FALSE))))) * wWeight[1]
    p$enr.score2 <- (unlist(viewSums(Views(x, flank(p, width=wSizes[2], start=TRUE))))+
                       unlist(viewSums(Views(x, flank(p, width=wSizes[2], start=FALSE))))) * wWeight[2]
    x <- p$weight * (p$enr.score1+p$enr.score2)
    c(mean=mean(x), sd=sd(x), len=length(p))
  })
  do.call(rbind, res)
}

runTFA <- function(pmoi, readlist, design){
  
  ptm <- proc.time()
  
  pmoi <- addMotifRelativeWeights(pmoi)
  se <- getMotifEnrScores(readlist, pmoi, BPPARAM=MulticoreParam(12))
  e <- log(assay(se)+0.1) # we ideally need to log-transform for limma, maybe double-check that the pseudo-count is appropriate
  fit <- eBayes(lmFit(e, design))
  TFA <- as.data.frame(topTable(fit, number = Inf))
  
  runtime <- proc.time()-ptm
  
  return(list(TFA, runtime))
}
