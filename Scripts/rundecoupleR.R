rundecoupleR <- function(counts,
                         genome,
                         motifs,
                         DAR,
                         decoupleR_modes=c("fgsea", 
                                           #"aucell", 
                                           "mlm", 
                                           "ulm", 
                                           "consensus", 
                                           "udt")){
  
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
  
  # Running all methods mentioned in decoupleR_modes
  
  # fgsea
  
  if ("fgsea" %in% decoupleR_modes){
    
    ptm <- proc.time()
    dcplr_gsea <- run_fgsea(DARmat, network, .source='source', .target='target', eps = 0, minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_gsea_raw <- list(dcplr_gsea, runtime)
    saveRDS(dcplr_gsea_raw, "./results/raw/dcplr_gsea_raw.rds")
    
    dcplr_gseadf <- as.data.frame(dcplr_gsea)
    dcplr_gseadf <- subset(dcplr_gseadf, statistic == "fgsea")
    dcplr_gseadf <- dcplr_gseadf[order(dcplr_gseadf$p_value),]
    dcplr_gseadf <- data.frame(row.names=dcplr_gseadf$source,
                               score = dcplr_gseadf$score,
                               padj = p.adjust(dcplr_gseadf$p_value, method="fdr"),
                               p = dcplr_gseadf$p_value)
    dcplr_gseadf$rank = seq_along(row.names(dcplr_gseadf))
    
    saveRDS(dcplr_gseadf, "./results/with_pvalues/dcplr_gsea.rds")
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
  #   saveRDS(dcplr_aucell_raw, "./results/raw/dcplr_aucell_raw.rds")
  # 
  #   dcplr_aucelldf <- as.data.frame(dcplr_aucell)
  #   dcplr_aucelldf <- dcplr_aucelldf[order(dcplr_aucelldf$score, decreasing = TRUE),]
  #   dcplr_aucelldf <- data.frame(row.names=dcplr_aucelldf$source,
  #                                score = dcplr_aucelldf$score)
  #   dcplr_aucelldf$rank = seq_along(row.names(dcplr_aucelldf))
  # 
  #   saveRDS(dcplr_aucelldf, "./results/scores_only/dcplr_aucell.rds")
  # }
  
  # udt
  
  if ("udt" %in% decoupleR_modes){
    
    ptm <- proc.time()
    dcplr_udt <- run_udt(DARmat, network, .source='source', .target='target', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_udt_raw <- list(dcplr_udt, runtime)
    saveRDS(dcplr_udt_raw, "./results/raw/dcplr_udt_raw.rds")
    
    dcplr_udtdf <- as.data.frame(dcplr_udt)
    dcplr_udtdf <- dcplr_udtdf[order(dcplr_udtdf$score, decreasing = TRUE),]
    dcplr_udtdf <- data.frame(row.names=dcplr_udtdf$source,
                                 score = dcplr_udtdf$score)
    dcplr_udtdf$rank = seq_along(row.names(dcplr_udtdf))
    
    saveRDS(dcplr_udtdf, "./results/scores_only/dcplr_udt.rds")
  }
  
  # mlm
  
  if ("mlm" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplr_mlm <- run_mlm(DARmat, network, .source='source', .target='target', .mor='mor', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_mlm_raw <- list(dcplr_mlm, runtime)
    saveRDS(dcplr_mlm_raw, "./results/raw/dcplr_mlm_raw.rds")
    
    dcplr_mlmdf <- as.data.frame(dcplr_mlm)
    dcplr_mlmdf <- dcplr_mlmdf[order(dcplr_mlmdf$p_value),]
    dcplr_mlmdf <- data.frame(row.names=dcplr_mlmdf$source,
                              score = dcplr_mlmdf$score,
                              padj = p.adjust(dcplr_mlmdf$p_value, method="fdr"),
                              p = dcplr_mlmdf$p_value)
    dcplr_mlmdf$rank = seq_along(row.names(dcplr_mlmdf))
    
    saveRDS(dcplr_mlmdf, "./results/with_pvalues/dcplr_mlm.rds")
  }
  
  # ulm
  
  if ("ulm" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplr_ulm <- run_ulm(DARmat, network, .source='source', .target='target', .mor='mor', minsize = 0)
    runtime <- proc.time()-ptm
    
    dcplr_ulm_raw <- list(dcplr_ulm, runtime)
    saveRDS(dcplr_ulm_raw, "./results/raw/dcplr_ulm_raw.rds")
    
    dcplr_ulmdf <- as.data.frame(dcplr_ulm)
    dcplr_ulmdf <- dcplr_ulmdf[order(dcplr_ulmdf$p_value),]
    dcplr_ulmdf <- data.frame(row.names=dcplr_ulmdf$source,
                              score = dcplr_ulmdf$score,
                              padj = p.adjust(dcplr_ulmdf$p_value, method="fdr"),
                              p=dcplr_ulmdf$p_value)
    dcplr_ulmdf$rank = seq_along(row.names(dcplr_ulmdf))
    
    saveRDS(dcplr_ulmdf, "./results/with_pvalues/dcplr_ulm.rds")
  } 
  
  # consensus
  
  if ("consensus" %in% decoupleR_modes){
    ptm <- proc.time()
    dcplR <- decouple(DARmat, network, .source='source', .target='target', minsize = 0)
    consensus <- run_consensus(dcplR)
    runtime <- proc.time()-ptm
    
    dcplr_consensus_raw <- list(consensus, runtime)
    saveRDS(dcplr_consensus_raw, "./results/raw/dcplr_consensus_raw.rds")
    
    dcplr_consensusdf <- as.data.frame(consensus)
    dcplr_consensusdf <- subset(dcplr_consensusdf, statistic == "consensus")
    dcplr_consensusdf <- subset(dcplr_consensusdf, condition == "1") # For some reason, decouple adds another condition called "V1" with different score and p-values.
    dcplr_consensusdf <- dcplr_consensusdf[order(dcplr_consensusdf$p_value),]
    dcplr_consensusdf <- data.frame(row.names = dcplr_consensusdf$source,
                              score = dcplr_consensusdf$score,
                              padj = p.adjust(dcplr_consensusdf$p_value, method="fdr"),
                              p = dcplr_consensusdf$p_value)
    dcplr_consensusdf$rank = seq_along(row.names(dcplr_consensusdf))
    
    saveRDS(dcplr_consensusdf, "./results/with_pvalues/dcplr_consensus.rds")
  } 
}