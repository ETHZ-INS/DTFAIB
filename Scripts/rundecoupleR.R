rundecoupleR <- function(counts,
                         genome,
                         motifs,
                         DAR,
                         decoupleR_modes="consensus"){
  
  # Generate required matrix of logFCs
  
  DARmat <- as.matrix(DAR$logFC)
  row.names(DARmat) <- rownames(DAR)
  
  # Generate required network
  
  matchingMotifs <- matchMotifs(motifs, counts, genome)
  row.names(matchingMotifs) <- as.character(granges(matchingMotifs))
  matchingMtx <- as(assay(matchingMotifs), "TsparseMatrix")
  network <- data.frame(source=factor(colnames(matchingMtx), colnames(matchingMtx))[matchingMtx@j+1L], 
                        target=factor(row.names(matchingMtx), row.names(matchingMtx))[matchingMtx@i+1L])
  network$mor <- 1L
  network$target <- sub(":\\+$", "", network$target)
  
  # remove colinearity
  network <- .removeColinearity(network)
  
  # Running all methods mentioned in decoupleR_modes
  
  # fgsea
  
  if ("fgsea" %in% decoupleR_modes){
    
    ptm <- proc.time()
    dcplr_gsea <- run_fgsea(DARmat, network, .source='source', .target='target', eps = 0, minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_gsea_raw <- list(dcplr_gsea, runtime)
    saveRDS(dcplr_gsea_raw, "./runATAC_results/raw/dRgsea_raw.rds")
    
    dcplr_gseadf <- as.data.frame(dcplr_gsea)
    dcplr_gseadf <- subset(dcplr_gseadf, statistic == "fgsea")
    dcplr_gseadf <- dcplr_gseadf[order(dcplr_gseadf$p_value, -abs(dcplr_gseadf$score)),]
    dcplr_gseadf <- data.frame(row.names=dcplr_gseadf$source,
                               score = dcplr_gseadf$score,
                               padj = p.adjust(dcplr_gseadf$p_value, method="fdr"),
                               p = dcplr_gseadf$p_value)
    dcplr_gseadf$rank = seq_along(row.names(dcplr_gseadf))
    
    saveRDS(dcplr_gseadf, "./runATAC_results/with_pvalues/dRgsea.rds")
  }
  
  # # AUCell
  # 
  # if ("aucell" %in% decoupleR_modes){
  # 
  #   ptm <- proc.time()
  #   dcplr_aucell <- decoupleR::run_aucell(DARmat, network, .source='source', .target='target', minsize=0, aucMaxRank = 3)
  #   runtime <- proc.time()-ptm
  # 
  #   dcplr_aucell_raw <- list(dcplr_aucell, runtime)
  #   saveRDS(dcplr_aucell_raw, "./results/raw/dRaucell_raw.rds")
  # 
  #   dcplr_aucelldf <- as.data.frame(dcplr_aucell)
  #   dcplr_aucelldf <- dcplr_aucelldf[order(dcplr_aucelldf$score, decreasing = TRUE),]
  #   dcplr_aucelldf <- data.frame(row.names=dcplr_aucelldf$source,
  #                                score = dcplr_aucelldf$score)
  #   dcplr_aucelldf$rank = seq_along(row.names(dcplr_aucelldf))
  # 
  #   saveRDS(dcplr_aucelldf, "./results/scores_only/dRaucell.rds")
  # }
  
  # udt
  
  if ("udt" %in% decoupleR_modes){
    
    ptm <- proc.time()
    dcplr_udt <- run_udt(DARmat, network, .source='source', .target='target', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_udt_raw <- list(dcplr_udt, runtime)
    saveRDS(dcplr_udt_raw, "./runATAC_results/raw/dRudt_raw.rds")
    
    dcplr_udtdf <- as.data.frame(dcplr_udt)
    dcplr_udtdf <- dcplr_udtdf[order(dcplr_udtdf$score, decreasing = TRUE),]
    dcplr_udtdf <- data.frame(row.names=dcplr_udtdf$source,
                                 score = dcplr_udtdf$score)
    dcplr_udtdf$rank = seq_along(row.names(dcplr_udtdf))
    
    saveRDS(dcplr_udtdf, "./runATAC_results/scores_only/dRudt.rds")
  }
  
  # mlm
  
  if ("mlm" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplr_mlm <- run_mlm(DARmat, network, .source='source', .target='target', .mor='mor', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_mlm_raw <- list(dcplr_mlm, runtime)
    saveRDS(dcplr_mlm_raw, "./runATAC_results/raw/dRmlm_raw.rds")
    
    dcplr_mlmdf <- as.data.frame(dcplr_mlm)
    dcplr_mlmdf <- dcplr_mlmdf[order(dcplr_mlmdf$p_value),]
    dcplr_mlmdf <- data.frame(row.names=dcplr_mlmdf$source,
                              score = dcplr_mlmdf$score,
                              padj = p.adjust(dcplr_mlmdf$p_value, method="fdr"),
                              p = dcplr_mlmdf$p_value)
    dcplr_mlmdf$rank = seq_along(row.names(dcplr_mlmdf))
    
    saveRDS(dcplr_mlmdf, "./runATAC_results/with_pvalues/dRmlm.rds")
  }
  
  # ulm
  
  if ("ulm" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplr_ulm <- run_ulm(DARmat, network, .source='source', .target='target', .mor='mor', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_ulm_raw <- list(dcplr_ulm, runtime)
    saveRDS(dcplr_ulm_raw, "./runATAC_results/raw/dRulm_raw.rds")
    
    dcplr_ulmdf <- as.data.frame(dcplr_ulm)
    dcplr_ulmdf <- dcplr_ulmdf[order(dcplr_ulmdf$p_value, -abs(dcplr_ulmdf$score)),]
    dcplr_ulmdf <- data.frame(row.names=dcplr_ulmdf$source,
                              score = dcplr_ulmdf$score,
                              padj = p.adjust(dcplr_ulmdf$p_value, method="fdr"),
                              p=dcplr_ulmdf$p_value)
    dcplr_ulmdf$rank = seq_along(row.names(dcplr_ulmdf))
    
    saveRDS(dcplr_ulmdf, "./runATAC_results/with_pvalues/dRulm.rds")
  } 
  
  # consensus
  
  if ("consensus" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplR <- decouple(DARmat, network, .source='source', .target='target', minsize = 0)
    consensus <- run_consensus(dcplR)
    runtime <- proc.time()-ptm
    
    dcplr_consensus_raw <- list(consensus, runtime)
    saveRDS(dcplr_consensus_raw, "./runATAC_results/raw/dRconsensus_raw.rds")
    
    dcplr_consensusdf <- as.data.frame(consensus)
    dcplr_consensusdf <- subset(dcplr_consensusdf, statistic == "consensus")
    dcplr_consensusdf <- subset(dcplr_consensusdf, condition == "1") # For some reason, decouple adds another condition called "V1" with different score and p-values.
    dcplr_consensusdf <- dcplr_consensusdf[order(dcplr_consensusdf$p_value),]
    dcplr_consensusdf <- data.frame(row.names = dcplr_consensusdf$source,
                              score = dcplr_consensusdf$score,
                              padj = p.adjust(dcplr_consensusdf$p_value, method="fdr"),
                              p = dcplr_consensusdf$p_value)
    dcplr_consensusdf$rank = seq_along(row.names(dcplr_consensusdf))
    
    saveRDS(dcplr_consensusdf, "./runATAC_results/with_pvalues/dRconsensus.rds")
  } 
}


.removeColinearity <- function(reg, minsize=5, truth=c(), cor.thres=0.96){
  # remove colinearity
  m <- reshape2::dcast(reg, target~source, value.var = "mor", fill = 0)
  row.names(m) <- m[,1]
  m <- as.matrix(m[,-1])
  m <- m[row.names(m),]
  m <- m[,colSums(m!=0)>=minsize]
  cc <- cor(m)
  cc[upper.tri(cc, diag=TRUE)] <- NA_real_
  removed <- vector("character")
  while(nrow(w <- which(cc>cor.thres, arr.ind=TRUE))>0){
    # if the true TF is colinear with something else, keep the true one
    isTrueReg <- intersect(row.names(w), truth)
    removed <- c(removed, setdiff(row.names(w),truth))
    if(length(isTrueReg)>0) removed <- c(removed, colnames(cc)[w[isTrueReg,2]])
    keep <- setdiff(row.names(cc),row.names(w))
    cc <- cc[keep,][,keep]
  }
  message(paste0(ncol(cc),"/",length(unique(reg$source))," regulons kept."))
  if(length(removed)>0)
    message("The following factors were removed due to collinearity with other factors:\n",
            paste(removed, collapse=", "))
  reg[reg$source %in% row.names(cc),]
}