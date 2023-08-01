suppressPackageStartupMessages({
  library(viper)
  library(SummarizedExperiment)
  library(Matrix)
  library(GenomicRanges)
  library(edgeR)
  library(limma)
  library(monaLisa)
  library(fgsea)
  library(stringr)
  library(chromVAR)
  library(motifmatchr)
  library(BiocParallel)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(Rsamtools)
  library(dplyr)
  library(MotifDb)
  library(universalmotif)
  library(memes)
  library(GenomicAlignments)
  library(data.table)
  library(decoupleR)
  library(AUCell)
  library(TFBSTools)
  library(ggplot2)
  library(epiwraps)
})

# Paths to the method wrappers:

source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runVIPER.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runGSEA.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runmonaLisa.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runchromVAR.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runSTOCKchromVAR.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/getpmoi.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runStabSel.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/rundecoupleR.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runregreg.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runulm.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runBaGFoot.R")

#' Title
#'
#' @param design vector of 1s and -1s to indicate the experimental design. E.g.: c(1,1,1,-1,-1,-1)
#' @param genome respective genome used for alignment
#' @param peakpath path to the ATAC-Seq peaks
#' @param regulons 
#' @param genesets 
#' @param readlist 
#' @param readtype 
#'
#' @return
#' @export
#'
#' @examples

runATAC <- function(methods=c("chromVAR", 
                              "monaLisa", 
                              "GSEA", 
                              "VIPER", 
                              "msVIPER", 
                              "diffTF", 
                              "StabSel", 
                              "decoupleR", 
                              "ulm",
                              "regreg",
                              "BaGFoot",
                              "MBA"),
                    decoupleR_modes=c("fgsea", 
                                      "mlm", 
                                      "ulm", 
                                      "consensus", 
                                      "udt"),
                    design,
                    genome,
                    peakpath,
                    pmoi,
                    spec=c("Hsapiens", "Mmusculus"),
                    readlist,
                    readtype=c("bam", "bed"),#[1+as.integer(grepl("bed$",readlist[1]))],
                    seqStyle=c("ensembl","UCSC"), # here the seqStyle refers to the style of the alignment files and not the peaks
                    diffTFpath,
                    rndSeed=1997,
                    counts=NULL,
                    paired_arg = TRUE)
{
  
  methods <- match.arg(methods, several.ok = TRUE)
  
  # Obtaining and re-sizing peaks
  set.seed(rndSeed)
  peakfile <- file(peakpath)
  peaks <- chromVAR::getPeaks(peakfile, 
                              sort_peaks = TRUE)
  peaks <- resize(peaks, 
                  width = 300, 
                  fix = "center")
  peaks <- keepStandardChromosomes(peaks, 
                                   pruning.mode = "coarse")
  
  # Obtaining read counts for the peaks and dividing them by their conditions
  
  seqlevelsStyle(peaks) <- seqStyle
  
  if ("chromVAR" %in% methods | "monaLisa" %in% methods | "GSEA" %in% methods | "VIPER" %in% methods | "msVIPER" %in% methods | "StabSel" %in% methods | "decoupleR" %in% methods | "ulm" %in% methods | "regreg" %in% methods) {
  
  if (!is.null(counts)) {
  # Skip countmatrix generation
  } else {  
    
  set.seed(rndSeed)
  counts <- chromVAR::getCounts(readlist,
                                peaks,
                                paired = paired_arg,
                                format = readtype)
  saveRDS(counts, "./runATAC_results/others/countmatrix.rds")
  }
    
  npos <- sum(design == 1)
  nneg <- sum(design == -1)
  counts_control <- counts[, colnames(counts)[1:npos]]
  counts_perturbed <- counts[, colnames(counts)[(npos+1):(npos+nneg)]]
  }
  
  seqlevelsStyle(pmoi) <- seqStyle
  
  pmoi2 <- pmoi # needed unmodified for BaGFoot
  
  # Compute regulons and GSEA genesets from pmoi object
  
  pmoi$motif_id <- factor(pmoi$motif_id)
  
    # divide motif scores by the max for that motif
  
  msmax <- max(splitAsList(pmoi$score, pmoi$motif_id))
  pmoi$score2 <- pmoi$score/msmax[pmoi$motif_id]
  
    # overlap peaks with motif instances
  
  o <- findOverlaps(peaks, pmoi, ignore.strand=TRUE)
  
    # build a matrix of scores (e.g. for alternative use with chromVAR)
  
  m <- sparseMatrix(from(o), as.integer(pmoi$motif_id[to(o)]), x=pmoi$score2[to(o)])
  row.names(m) <- as.character(granges(peaks))
  colnames(m) <- levels(pmoi$motif_id)
  
    # build a regulon (for use with viper)
    # here we need to name the peaks (I use the coordinates) and make sure the names
    #  corresponding to those of the differential accessibility signature later on
  
  ll <- split( data.frame(peak=as.character(granges(peaks))[from(o)], score=pmoi$score2[to(o)]),
               pmoi$motif_id[to(o)] )
  
  regulons <- lapply(ll, FUN=function(x){
    # count a peak only once if it has two instances of the same motif
    x <- x[order(-x$score),]
    x <- x[!duplicated(x$score),]
    # reformat to regulon
    list( tfmode=setNames(rep(1L,nrow(x)), x$peak),
          likelihood=setNames(x$score, x$peak) )
  })
  
  # Changing the names in the regulons to the conventional ones
  
  if (spec=="Hsapiens") {
    motifnames <- fread("/mnt/plger/fgerbaldo/DTFAIB/Scripts/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv")
  } else if (spec=="Mmusculus") {
    motifnames <- fread("/mnt/plger/fgerbaldo/DTFAIB/Scripts/HOCOMOCOv11_core_annotation_MOUSE_mono.tsv")
  }
  
  # These will be turned into the motifs used for chromVAR and monaLisa later but first, they serve as donor of "wrong" TF names
  
  motifs <- getNonRedundantMotifs(format = "PWMatrix", 
                                  species = spec)
  
  # Listing TF names we don't want (wrong)
  
  motiflist <- c()
  
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        motiflist$wrong <- append(motiflist$wrong,motifs[[i]]@name)
    }
  }
  
  # Adding a list entry with the TF names we want (correct)
  
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        motiflist$correct <- append(motiflist$correct,motifnames$`Transcription factor`[[j]])
    }
  }
  
  # Replacing the undesired names in the regulons with the desired ones
  
  for (i in seq_along(regulons)){
    for (j in seq_along(motiflist$wrong)){
      if (names(regulons)[[i]] == motiflist$wrong[[j]])
        names(regulons)[[i]] <- motiflist$correct[[j]]
    }
  }
  
  saveRDS(regulons, "./runATAC_results/others/regulons.rds")
  
  # get normal lists of peaks per motif (for GSEA)
  
  genesets <- split( as.character(granges(peaks))[from(o)], pmoi$motif_id[to(o)] )
  
  # Replacing the undesired names in the genesets with the desired ones
  
  for (i in seq_along(genesets)){
    for (j in seq_along(motiflist$wrong)){
      if (names(genesets)[[i]] == motiflist$wrong[[j]])
        names(genesets)[[i]] <- motiflist$correct[[j]]
    }
  }
  
  saveRDS(genesets, "./runATAC_results/others/genesets.rds")
  
  # Here, the names in the motifs for chromVAR and monaLisa are changed to the correct ones, as well
    
  for (i in seq_along(motifs)){
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@ID == motifnames$Model[[j]])
        names(motifs)[[i]] <- motifnames$`Transcription factor`[[j]]
    }
  }
  
  BANP_motif <- readRDS("/mnt/plger/fgerbaldo/DTFAIB/Scripts/BANP.PFMatrix.rds")
  BANP_motifPW <- TFBSTools::toPWM(BANP_motif) 
  motifs$BANP <- BANP_motifPW
  
  saveRDS(motifs, "./runATAC_results/others/motifs.rds")
  
  seqlevelsStyle(genome) <- seqStyle
  
  # Compute differentially accessible regions required to run monaLisa, StabSel, fGSEA, VIPER, and msVIPER 
  
  if ("monaLisa" %in% methods | "StabSel" %in% methods | "GSEA" %in% methods | "VIPER" %in% methods | "msVIPER" %in% methods | "decoupleR" %in% methods | "ulm" %in% methods | "regreg" %in% methods) {
    
    set.seed(rndSeed)
    DAR <- dATestedgeR(counts_control, 
                       counts_perturbed)
    
    # Generate required matrix of logFCs
    
    DARmat <- as.numeric(DAR$logFC)
    names(DARmat) <- rownames(DAR)
  }
  
  if ("ulm" %in% methods | "regreg" %in% methods) {
    set.seed(rndSeed)
    matchingMotifs <- matchMotifs(motifs, counts, genome)
    row.names(matchingMotifs) <- as.character(granges(matchingMotifs))
    matchMtx <- as(assay(matchingMotifs), "dgCMatrix")
    
    rownames(matchMtx) <- gsub(":\\+$", "", rownames(matchMtx))
  }
  
  # Start running methods
  
  readouts <- list()
  
  # Run chromVAR without normalization
  
  if ("chromVAR" %in% methods){
    set.seed(rndSeed)
    CV <- runSTOCKchromVAR(counts, 
                      genome, 
                      motifs, 
                      design)
    
    saveRDS(CV, "./runATAC_results/raw/CV_raw.rds")
    
    # Select desired information from chromVAR readout
    
    CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
    CVdf <- data.frame(row.names = rownames(CVsel),
                       logFC = CVsel[, "logFC"], 
                       padj = CVsel[, "adj.P.Val"], 
                       p = CVsel[, "P.Value"])
    CVdf$rank = seq_along(row.names(CVdf))
    
    saveRDS(CVdf, "./runATAC_results/with_pvalues/CV.rds")
    readouts$CV <- CVdf
  }
  
  # Run chromVAR with normalization
  
  if ("chromVAR" %in% methods){
  set.seed(rndSeed)
  CV <- runchromVAR(counts, 
                    genome, 
                    motifs, 
                    design)
  
  saveRDS(CV, "./runATAC_results/raw/CVnorm_raw.rds")
  
  # Select desired information from chromVAR readout
  
  CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
  CVdf <- data.frame(row.names = rownames(CVsel),
                     logFC = CVsel[, "logFC"], 
                     padj = CVsel[, "adj.P.Val"], 
                     p = CVsel[, "P.Value"])
  CVdf$rank = seq_along(row.names(CVdf))
  
  saveRDS(CVdf, "./runATAC_results/with_pvalues/CVnorm.rds")
  readouts$CVnorm <- CVdf
  }
  
  if ("decoupleR" %in% methods){
  
  set.seed(rndSeed)
  rundecoupleR(counts,
               genome,
               motifs,
               DAR,
               decoupleR_modes = decoupleR_modes)
  }
  
  # Run monaLisa 
  
  if ("monaLisa" %in% methods){
    set.seed(rndSeed)
  ML <- runmonaLisa(DAR, 
                    motifs, 
                    peaks, 
                    genome)
  
  saveRDS(ML, "./runATAC_results/raw/ML_raw.rds")
  
    # Calculate bin-level p-values based on Simes method with code from https://github.com/markrobinsonuzh/DAMEfinder/blob/master/R/simes_pval.R
  
  simes <- function(pval){ 
    min((length(pval)*pval[order(pval)])/seq(from = 1, to = length(pval), by = 1))
  }
  
  ML <- ML[[1]]
  ML <- ML[, colData(ML)$bin.nochange == FALSE]
  MLp <- assays(ML)$negLog10P
  MLp <- 10^(MLp*-1)
  
  MLsimes <- apply(MLp, 1, simes)
  
  # Correct the p-values using FDR correction 
  
  MLdf <- data.frame(p = sort(MLsimes))
  MLdf$rank <- seq_along(row.names(MLdf))
  MLpadj <- p.adjust(MLdf$p, method="fdr")
  MLdf <- cbind(padj = p.adjust(MLdf$p, method="fdr"), MLdf)
  
  saveRDS(MLdf, "./runATAC_results/with_pvalues/ML.rds")
  readouts$ML <- MLdf
  
  # calculate correlation across bins
  # cors <- cor(t(assays(ML[[1]])$log2enr), seq_len(ncol(ML[[1]])), method="spearman")[,1]
  # names(cors) <- row.names(ML[[1]])
  # MLdf$binSpearman <- cors[row.names(MLdf),]
  # MLdf <- MLdf[order(abs(MLdf$binSpearman)*-log10(MLdf$p), decreasing=TRUE),]
  # saveRDS(MLdf, "./runATAC_results/with_pvalues/MLsp.rds")
  # readouts$MLsp <- MLdf
  
  }
  
  # Run StabSel
  
  if ("StabSel" %in% methods){
    set.seed(rndSeed)
  MLStabSel <- runStabSel(DAR, 
                      motifs, 
                      peaks, 
                      genome)
  
  saveRDS(MLStabSel, "./runATAC_results/raw/MLStabSel_raw.rds")
  
  MLStabSeldf <- data.frame(TF = colnames(MLStabSel[[1]]),
                         sel_Prob = MLStabSel[[1]]$selProb)
  MLStabSeldf <- MLStabSeldf[order(MLStabSeldf$sel_Prob, decreasing = TRUE),]         
  MLStabSeldf <- data.frame(row.names = MLStabSeldf$TF,
                         sel_Prob = MLStabSeldf$sel_Prob)
  MLStabSeldf <- subset(MLStabSeldf, !(row.names(MLStabSeldf) %in% c("fracGC", "oeCpG")))
  MLStabSeldf$rank <- seq_along(rownames(MLStabSeldf))
  
  saveRDS(MLStabSeldf, "./runATAC_results/scores_only/MLStabSel.rds")
  readouts$MLStabSel <- MLStabSeldf
  }
  
  # Run fGSEA
  
  if ("GSEA" %in% methods){
  
  set.seed(rndSeed)
  GSEA <- rungsea(DAR, 
                  genesets, 
                  peaks)
  
  saveRDS(GSEA, "./runATAC_results/raw/GSEA_raw.rds")
  
  # Select desired information from fgsea readout
  
  GSEA <- GSEA[[1]][order(GSEA[[1]]$pval, -abs(GSEA[[1]]$NES)),]
  GSEAdf <- data.frame(row.names = GSEA$pathway,
                       NES = GSEA$NES,
                       padj = GSEA$padj,
                       p = GSEA$pval)
  
  GSEAdf$rank <- seq_along(row.names(GSEAdf))
  
  saveRDS(GSEAdf, "./runATAC_results/with_pvalues/GSEA.rds")
  readouts$GSEA <- GSEAdf
  }
  
  # Run VIPER
  if ("VIPER" %in% methods){
  
  set.seed(rndSeed)
  VIPER <- runviper(counts_control, 
                    counts_perturbed, 
                    DAR, 
                    peaks, 
                    regulons, 
                    mode = "viper", 
                    method = "ttest")
  
  saveRDS(VIPER, "./runATAC_results/raw/VIPER_raw.rds")
  
  VIPERtTest <- rowTtest(VIPER[[1]][,1:npos], VIPER[[1]][,(npos+1):(npos+nneg)])
  
  # Select desired information from VIPER readout
  
  VIPERdf <- data.frame(VIPERtTest)
  VIPERdf <- VIPERdf[order(VIPERdf$p.value),]
  VIPERdf$padj <- p.adjust(VIPERdf$p.value, method="fdr")
  VIPERdf <- data.frame(row.names = rownames(VIPERdf),
                        statistic = VIPERdf$statistic,
                        padj = VIPERdf$padj,
                        p = VIPERdf$p.value)
  VIPERdf$rank <- seq_along(row.names(VIPERdf))
  
  saveRDS(VIPERdf, "./runATAC_results/with_pvalues/VIPER.rds")
  readouts$VIPER <- VIPERdf
  }
  
  # Run msVIPER
  
  if ("msVIPER" %in% methods){
  set.seed(rndSeed)
  msVIPER <- runviper(counts_control, 
                      counts_perturbed,
                      DAR,
                      peaks,
                      regulons,
                      mode = "msviper",
                      method = "ttest")
  
  saveRDS(msVIPER, "./runATAC_results/raw/msVIPER_raw.rds")
  
  # Select desired information from msVIPER readout
  
  msVIPERdf <- data.frame(NES = msVIPER[[1]]$es$nes,
                          padj = p.adjust(msVIPER[[1]]$es$p.value, method="fdr"),
                          p = msVIPER[[1]]$es$p.value)
  
  msVIPERdf <- msVIPERdf[order(msVIPERdf$p),]
  msVIPERdf$rank <- seq_along(row.names(msVIPERdf))
  
  saveRDS(msVIPERdf, "./runATAC_results/with_pvalues/msVIPER.rds")
  readouts$msVIPER <- msVIPERdf
  }
  
  # Run ulm
  
  if ("ulm" %in% methods){
    set.seed(rndSeed)
    ulm <- runulm(DARmat = DARmat,
                  matchMtx = matchMtx)
    
    saveRDS(ulm, "./runATAC_results/raw/ulm_raw.rds")
    
    ulmdf <- data.frame(row.names = rownames(ulm[[1]]),
                        score = ulm[[1]]$score,
                        padj = ulm[[1]]$padj,
                        p = ulm[[1]]$p)
    ulmdf <- ulmdf[order(ulmdf$p, -abs(ulmdf$score)),]
    ulmdf$rank <- seq_along(row.names(ulmdf))
    
    
    saveRDS(ulmdf, "./runATAC_results/with_pvalues/ulm.rds")
    readouts$ulm <- ulmdf
  }
  
  # Run regreg
  
  if ("regreg" %in% methods){
    set.seed(rndSeed)
    regreg <- runregreg(DARmat = DARmat,
                        matchMtx = matchMtx)
    
    saveRDS(regreg, "./runATAC_results/raw/regreg_raw.rds")
    
    regregdf <- data.frame(row.names = rownames(regreg[[1]]),
                           score = regreg[[1]]$score,
                           padj = regreg[[1]]$FDR,
                           p = regreg[[1]]$p_value)
    regregdf$rank <- seq_along(row.names(regregdf))
    
    saveRDS(regregdf, "./runATAC_results/with_pvalues/regreg.rds")
    readouts$regreg <- regregdf
  }
  
  # # Run PLG's TFA
  # 
  # if ("TFA" %in% methods){
  # set.seed(rndSeed)
  # TFA <- runTFA(pmoi,
  #               readlist,
  #               design)
  # 
  # saveRDS(TFA, "./runATAC_results/raw/TFA_raw.rds")
  # 
  # # Select desired information from TFA readout
  #     
  # TFAdf <- data.frame(row.names = rownames(TFA[[1]]),
  #                     logFC = TFA[[1]]$logFC,
  #                     padj = TFA[[1]]$adj.P.Val,
  #                     p = TFA[[1]]$P.Value)
  # 
  # TFAdf$rank = seq_along(row.names(TFAdf))
  # 
  # saveRDS(TFAdf, "./runATAC_results/with_pvalues/TFA.rds")
  # } 
  
  # run BaGFootLike method
  
  if ("BaGFoot" %in% methods){
    set.seed(rndSeed)
    BaGFoot <- runBaGFoot(readlist,
                          pmoi2,
                          paired_arg)
    saveRDS(BaGFoot, "./runATAC_results/raw/BaGFootLike_raw.rds")
    
    BaGFoot <- BaGFoot[[1]]
    BFdf <- data.frame(row.names = rownames(BaGFoot),
                       padj=BaGFoot$FDR,
                       p=BaGFoot$combined.P)
    BFdf <- BFdf[order(BFdf$p),]
    BFdf$rank = seq_along(row.names(BFdf))
    saveRDS(BFdf,"./runATAC_results/with_pvalues/BaGFootLike.rds")
    readouts$BaGFoot <- BFdf
  }
  
  # run Model-based analysis
  
  if ("MBA" %in% methods){
    set.seed(rndSeed)
    MBA <- runMBA(readlist,
                  pmoi2,
                  paired_arg)
    saveRDS(MBA, "./runATAC_results/raw/MBA_raw.rds")
    MBA <- MBA[[1]]
    MBAdf <- data.frame(row.names=rownames(MBA),
                        padj=MBA$adj.P.Val,
                        p=MBA$P.Value)
    MBAdf <- MBAdf[order(MBAdf$p),]
    MBAdf$rank =seq_along(row.names(MBAdf))
    saveRDS(MBAdf, "runATAC_results/with_pvalues/MBA.rds")
    readouts$MBA <- MBAdf
  }
  
  # Process pre-generated diffTF output
  
  if ("diffTF" %in% methods){
    set.seed(rndSeed)
  diffTFdata <- fread(diffTFpath)
  diffTFdf <- data.frame(row.names = diffTFdata$TF,
                           WMD = diffTFdata$weighted_meanDifference,
                           padj = diffTFdata$pvalueAdj,
                           p = diffTFdata$pvalue)
  diffTFdf <- diffTFdf[order(diffTFdf$p, -abs(diffTFdf$WMD)),]
  diffTFdf$rank <- seq_along(row.names(diffTFdf))
  saveRDS(diffTFdf, "./runATAC_results/with_pvalues/diffTF.rds")
  readouts$diffTF <- diffTFdf
  }
  saveRDS(readouts, "./runATAC_results/others/readouts.rds")
  return(readouts=readouts)
}