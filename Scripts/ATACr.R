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
  library(metap)
  library(TFBSTools)
  library(ggplot2)
  library(epiwraps)
})

# Paths to the method wrappers:

source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runVIPER.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runGSEA.R") #still source this script for the future, eventhough the method is commented out due to NA values disrupting aggregation
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runmonaLisa.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runchromVAR.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runSTOCKchromVAR.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/getpmoi.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runulm.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/runregreg.R")
source("/mnt/plger/fgerbaldo/DTFAIB/Scripts/npPvalAgg.R")

ATACr <- function(methods=c("chromVAR", 
                               "monaLisa",
                               "msVIPER", 
                               "ulm", 
                               "regreg"),
                     design,
                     genome,
                     peakpath,
                     pmoi,
                     spec=c("Hsapiens", "Mmusculus"),
                     readlist,
                     readtype=c("bam", "bed"),#[1+as.integer(grepl("bed$",readlist[1]))],
                     seqStyle=c("ensembl","UCSC"), # here the seqStyle refers to the style of the alignment files and not the peaks
                    aggregation = FALSE, # pvalue aggregation or not?
                     rndSeed=1997,
                  counts=NULL, # if a countmatrix has been generated before, it can be provided to skip this step
                  paired_arg = TRUE) # paired reads or not?
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
  
  if (!is.null(counts)) {
    # Skip countmatrix generation
  } else {
  
  set.seed(rndSeed)
  counts <- chromVAR::getCounts(readlist,
                                peaks,
                                paired = paired_arg,
                                format = readtype)
  saveRDS(counts, "./ATACr_results/others/countmatrix.rds")
  }
  
  npos <- sum(design == 1)
  nneg <- sum(design == -1)
  counts_control <- counts[, colnames(counts)[1:npos]]
  counts_perturbed <- counts[, colnames(counts)[(npos+1):(npos+nneg)]]
  
  # Compute regulons and GSEA genesets from pmoi object
  
  seqlevelsStyle(pmoi) <- seqStyle
  
  pmoi2 <- pmoi # needed unmodified for BaGFoot
  
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
  
  saveRDS(regulons, "./ATACr_results/others/regulons.rds")
  
  # get normal lists of peaks per motif (for GSEA)
  
  genesets <- split( as.character(granges(peaks))[from(o)], pmoi$motif_id[to(o)] )
  
  # Replacing the undesired names in the genesets with the desired ones
  
  for (i in seq_along(genesets)){
    for (j in seq_along(motiflist$wrong)){
      if (names(genesets)[[i]] == motiflist$wrong[[j]])
        names(genesets)[[i]] <- motiflist$correct[[j]]
    }
  }
  
  saveRDS(genesets, "./ATACr_results/others/genesets.rds")
  
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
  
  saveRDS(motifs, "./ATACr_results/others/motifs.rds")
  
  seqlevelsStyle(genome) <- seqStyle
  
  # Compute differentially accessible regions required to run monaLisa, monaLasso, fGSEA, VIPER, and msVIPER 
  
  if ("monaLisa" %in% methods | "GSEA" %in% methods | "msVIPER" %in% methods | "ulm" %in% methods | "regreg" %in% methods) {
    
    set.seed(rndSeed)
    DAR <- dATestedgeR(counts_control, 
                       counts_perturbed)
  
  # Generate required matrix of logFCs
  
  DARmat <- as.numeric(DAR$logFC)
  names(DARmat) <- rownames(DAR)
  }
  
  # Generate required network for ulm and regreg
  
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
    
    saveRDS(CV, "./ATACr_results/raw/CV_raw.rds")
    
    # Select desired information from chromVAR readout
    
    CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
    CVdf <- data.frame(row.names = rownames(CVsel),
                       logFC = CVsel[, "logFC"], 
                       padj = CVsel[, "adj.P.Val"], 
                       p = CVsel[, "P.Value"])
    CVdf$rank = seq_along(row.names(CVdf))
    
    saveRDS(CVdf, "./ATACr_results/with_pvalues/CV.rds")
  }
  
  # Run chromVAR with normalization
  
  if ("chromVAR" %in% methods){
    set.seed(rndSeed)
    CV <- runchromVAR(counts, 
                      genome, 
                      motifs, 
                      design)
    
    saveRDS(CV, "./ATACr_results/raw/CVnorm_raw.rds")
    
    # Select desired information from chromVAR readout
    
    CVsel <- CV[[1]][order(CV[[1]]$adj.P.Val),]
    CVdf <- data.frame(row.names = rownames(CVsel),
                       logFC = CVsel[, "logFC"], 
                       padj = CVsel[, "adj.P.Val"], 
                       p = CVsel[, "P.Value"])
    CVdf$rank = seq_along(row.names(CVdf))
    
    saveRDS(CVdf, "./ATACr_results/with_pvalues/CVnorm.rds")
    readouts$CVnorm <- CVdf
  }
  
  # Run monaLisa 
  
  if ("monaLisa" %in% methods){
    set.seed(rndSeed)
    ML <- runmonaLisa(DAR, 
                      motifs, 
                      peaks, 
                      genome)
    
    saveRDS(ML, "./ATACr_results/raw/ML_raw.rds")
    
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
    
    saveRDS(MLdf, "./ATACr_results/with_pvalues/ML.rds")
    readouts$ML <- MLdf
    
    # calculate correlation across bins
    cors <- cor(t(assays(ML[[1]])$log2enr), seq_len(ncol(ML[[1]])), method="spearman")[,1]
    names(cors) <- row.names(ML[[1]])
    MLdf$binSpearman <- cors[row.names(MLdf),]
    MLdf <- MLdf[order(abs(MLdf$binSpearman)*-log10(MLdf$p), decreasing=TRUE),]
    saveRDS(MLdf, "./ATACr_results/with_pvalues/MLsp.rds")
    readouts$MLsp <- MLdf
  }
  
  # # Run fGSEA
  # 
  # if ("GSEA" %in% methods){
  #   
  #   set.seed(rndSeed)
  #   GSEA <- rungsea(DAR, 
  #                   genesets, 
  #                   peaks)
  #   
  #   saveRDS(GSEA, "./ATACr_results/raw/GSEA_raw.rds")
  #   
  #   # Select desired information from fgsea readout
  #   
  #   GSEA <- GSEA[[1]][order(GSEA[[1]]$pval, -abs(GSEA[[1]]$NES)),]
  #   GSEAdf <- data.frame(row.names = GSEA$pathway,
  #                        NES = GSEA$NES,
  #                        padj = GSEA$padj,
  #                        p = GSEA$pval)
  #   GSEAdf$rank <- seq_along(row.names(GSEAdf))
  #   
  #   saveRDS(GSEAdf, "./ATACr_results/with_pvalues/GSEA.rds")
  #   readouts$GSEA <- GSEAdf
  # }
  
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
    
    saveRDS(msVIPER, "./ATACr_results/raw/msVIPER_raw.rds")
    
    # Select desired information from msVIPER readout
    
    msVIPERdf <- data.frame(NES = msVIPER[[1]]$es$nes,
                            padj = p.adjust(msVIPER[[1]]$es$p.value, method="fdr"),
                            p = msVIPER[[1]]$es$p.value)
    
    msVIPERdf <- msVIPERdf[order(msVIPERdf$p),]
    msVIPERdf$rank <- seq_along(row.names(msVIPERdf))
    
    saveRDS(msVIPERdf, "./ATACr_results/with_pvalues/msVIPER.rds")
    readouts$msVIPER <- msVIPERdf
  }
  
  # Run ulm
  
  if ("ulm" %in% methods){
    set.seed(rndSeed)
    ulm <- runulm(DARmat = DARmat,
                  matchMtx = matchMtx)
    
    saveRDS(ulm, "./ATACr_results/raw/ulm_raw.rds")
    
    ulmdf <- data.frame(row.names = rownames(ulm[[1]]),
                        score = ulm[[1]]$score,
                        padj = ulm[[1]]$padj,
                        p = ulm[[1]]$p)
    ulmdf <- ulmdf[order(ulmdf$p, -abs(ulmdf$score)),]
    ulmdf$rank <- seq_along(row.names(ulmdf))
    
    
    saveRDS(ulmdf, "./ATACr_results/with_pvalues/ulm.rds")
    readouts$ulm <- ulmdf
  }
  
  # Run regreg
  
  if ("regreg" %in% methods){
    set.seed(rndSeed)
    regreg <- runregreg(DARmat = DARmat,
                        matchMtx = matchMtx)
    
    saveRDS(regreg, "./ATACr_results/raw/regreg_raw.rds")
    
    regregdf <- data.frame(row.names = rownames(regreg[[1]]),
                           score = regreg[[1]]$score,
                           padj = regreg[[1]]$FDR,
                           p = regreg[[1]]$p_value)
    regregdf$rank <- seq_along(row.names(regregdf))
    
    saveRDS(regregdf, "./ATACr_results/with_pvalues/regreg.rds")
    readouts$regreg <- regregdf
  }
  
  saveRDS(readouts, "./ATACr_results/others/readouts.rds")
  
  # p value aggregation (At this point only reasonable if all 5 methods are executed since the aggregation is always applied to the first 5 columns)
  
  if (aggregation){
    set.seed(rndSeed)
  # Since come methods drop some TFs entirely, we first need to determine a set of TFs that is contained in all readouts
  
  common_names <- names(motifs)
  for (df in readouts){
    common_names <- intersect(common_names, rownames(df))
  }
  
  # Generate a df with the adjusted pvalues 
  
  common_df <- data.frame(sample = common_names)
  for (i in seq_along(readouts)) {
    df <- readouts[[i]]
    df_common <- df[common_names, "p", drop = FALSE]
    colnames(df_common) <- names(readouts)[i]
    common_df <- cbind(common_df, df_common)
  }
  
  # get rid of the "sample" column again since now the TF names are in the rownames & set NA values from GSEA to 1
  
  # common_df[is.na(common_df)] <- 1 Bring this back when we decide to use GSEA
  
  common_df <- as.matrix(common_df[, -1])
  
  common_df[common_df< 1e-300] <- 1e-300
  common_df <- as.data.frame(common_df)
  
  # Pval Aggregation with Fisher's method
  
  common_df$fisher <- apply(common_df[, 1:5], 1, function(x) sumlog(x)$p)
  fisherdf <- data.frame(row.names = rownames(common_df),
                         padj = p.adjust(common_df$fisher, method="fdr"),
                         p = common_df$fisher)
  fisherdf <- fisherdf[order(fisherdf$p),]
  fisherdf$rank <- seq_along(row.names(fisherdf))
  saveRDS(fisherdf, "./ATACr_results/with_pvalues/fisher.rds")
  common_df$fisher <- p.adjust(common_df$fisher, method="fdr")
  
  # Pval Aggregation with Stouffer's method
  
  common_df$stouffer <- apply(common_df[, 1:5], 1, function(x) sumz(x)$p)
  stoufferdf <- data.frame(row.names = rownames(common_df),
                         padj = p.adjust(common_df$stouffer, method="fdr"),
                         p = common_df$stouffer)
  stoufferdf <- stoufferdf[order(stoufferdf$p),]
  stoufferdf$rank <- seq_along(row.names(stoufferdf))
  saveRDS(stoufferdf, "./ATACr_results/with_pvalues/stouffer.rds")
  common_df$stouffer <- p.adjust(common_df$stouffer, method="fdr")
  
  # Pval Aggregation with Simes' method
  
  common_df$simes <- apply(common_df[, 1:5], 1, simes)
  simesdf <- data.frame(row.names = rownames(common_df),
                           padj = p.adjust(common_df$simes, method="fdr"),
                           p = common_df$simes)
  simesdf <- simesdf[order(simesdf$p),]
  simesdf$rank <- seq_along(row.names(simesdf))
  saveRDS(simesdf, "./ATACr_results/with_pvalues/simes.rds")
  common_df$simes <- p.adjust(common_df$simes, method="fdr")
  
  # Pval Aggregation with a method for non-parametric aggregation
  
  npAggMat <- as.matrix(common_df[, 1:5])
  npPAggdf <- npPvalAgg(npAggMat, trans=c("sqrtrank"))
  saveRDS(npPAggdf, "./ATACr_results/with_pvalues/npPvalAgg.rds")
  common_df$npPvalAgg <- npPAggdf[rownames(common_df), "padj"]
  
  saveRDS(common_df, "./ATACr_results/summary.rds")
  }
}