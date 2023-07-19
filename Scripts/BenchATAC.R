library(SummarizedExperiment)
library(viper)
library(ggplot2)
library(tidyr)
source("/mnt/plger/fgerbaldo/BenchmarkTFactivity/BMScripts/getInteractors.R")

BenchATAC <- function(rawpaths, 
                      respaths,
                      rank_only,
                      TrueTF, 
                      motifpath, 
                      species=c("9606","10090")){ # species 9606 is for human and 10090 for mouse
  
  # make a list of all the readouts
  res <- list()
  for (path in respaths) {
    method_name <- tools::file_path_sans_ext(basename(path))
    res[[method_name]] <- readRDS(path)
  }
  
  # make a list of all the raw readouts which contain the run time informations
  raw.res <- list()
  for (rawpath in rawpaths) {
    method_name <- tools::file_path_sans_ext(basename(rawpath))
    raw.res[[method_name]] <- readRDS(rawpath)
  }
  
  # make a list of all the readouts from methods that provide no p-values but only scores and ranks
  
  if (is.null(rank_only)) {
    # execute code when no path is provided
    rank_res <- list()
  } else {
    # execute code when paths are provided
    rank_res <- list()
    for (path in rank_only) {
      method_name <- tools::file_path_sans_ext(basename(path))
      rank_res[[method_name]] <- readRDS(path)
    }
  }
  
  # rank_res <- list()
  # for (path in rank_only) {
  #   method_name <- tools::file_path_sans_ext(basename(path))
  #   rank_res[[method_name]] <- readRDS(path)
  # }
  
  # fetch the run time information
  runtime <- lapply(raw.res, function(x){
    print(x[[2]])
  })
  
  # combine methods with only ranks and those with both ranks and p-values
  rank_list <- c(res, rank_res)
  
  # fetch where the true transcription factor ranks according to its p-value
  true_rank <- lapply(rank_list, TrueTF, FUN=function(x, TrueTF){
    x$rank[which(row.names(x) %in% TrueTF)]
  })
  
  # put it in a format where the TF names are included:
  true_rank_table <- lapply(rank_list, TrueTF, FUN=function(x, TrueTF){
    cbind(rank = x$rank[which(row.names(x) %in% TrueTF)],
          name = rownames(x)[which(row.names(x) %in% TrueTF)])
  })
  
  # deal with integer(0) values
  true_rank <- lapply(true_rank, function(x) ifelse(length(x) == 0, NA, x))
  
  # make df of the true_ranks for plotting
  rank_df <- as.data.frame(true_rank)
  
  # reshaping the df because ggplot2 is not able to understand the data otherwise
  ranks <- gather(rank_df, method, rank)
  
  # Create the scatter plot
  rankplot <-  ggplot(ranks, aes(x = method, y = rank)) +
    scale_y_sqrt() +
    geom_point() +
    xlab("Method") +
    ylab("Rank")
  
  # fetch the adjusted p-value of the true transcription factor
  true_padj <- lapply(res, TrueTF, FUN=function(x, TrueTF){
    x[which(row.names(x) %in% TrueTF),"padj"]
  })
  
  # make df of the true_padjs to plot
  padj_df <- as.data.frame(true_padj)
  
  # reshaping the df because ggplot2 is not able to comprehend the data otherwise
  padj <- gather(padj_df, method, padj)
  padj$NegLogPadj <- -log(padj$padj)
  
  # plot the -log10(padj)
  padjplot <-  ggplot(padj, aes(x = method, y = NegLogPadj)) +
    geom_point() +
    xlab("Method") +
    ylab("-log(padj)")
  
  # fetch Interactors that were also part of the transcription factor list provided to each method
  Interactors <- unique(unlist(lapply(TrueTF, function(x) getInteractors(x, species))))
  motifnames <- readRDS(motifpath)
  motifnames <- names(motifnames)
  cofactors <- intersect(motifnames, Interactors)
  
  # fetch the number of significant transcription factors that are not the true TF
  `%!in%` = Negate(`%in%`)
  NotTrueTF <- lapply(res, TrueTF, FUN = function(x, TrueTF){
    sum(ifelse(x$padj<0.05 & row.names(x) %!in% TrueTF, 1, 0), na.rm = TRUE)
  })
  
  # fetch the number of significant transcription factors that are neither the trueTF nor one of its cofactors
  FD <- lapply(res, TrueTF, cofactors, FUN = function(x, TrueTF, cofactors){
    sum(ifelse(x$padj<0.05 & row.names(x) %!in% TrueTF & row.names(x) %!in% cofactors, 1, 0), na.rm = TRUE)
  })
  
  # Add all elements of trueTF to cofactors for AUC computation
  cofactors <- union(cofactors, TrueTF)
  
  res <- rank_list
  
  # add the cumulative sum of cofactors column to the results data frame
  res <- lapply(res, cofactors, FUN=function(x, cofactors){
    TruePositives <-  cumsum(as.integer(row.names(x) %in% cofactors))
    cbind(x, TruePositives)
  })
  
  # add the AUC_score column to the results data frame to compute the AUC
  res <- lapply(res, FUN=function(x){
    AUC_score <- x$TruePositives / x$rank
    cbind(x, AUC_score)
  })
  
  # This is a faster way to compute the AUC which however comes without being able to make a plot out of it. Thus, I leave both alternatives in for now.
  # res <- lapply(res, FUN=function(x){
  #   cofactor_logic <- rownames(x) %in% cofactors
  #   cbind(x, cofactor_logic)
  # })
  # AUC <- lapply(res, FUN=function(x){
  #   print(paste("AUC:", round(sum(cumsum(x$cofactor_logic) / seq_along(x$cofactor_logic)) / nrow(x), 2)))
  # })
  
  # compute the AUC based on the top 100 ranked TFs
  AUC <- lapply(res, FUN=function(x){
    print(round(mean(x$AUC_score[1:100]),2))
  })
  
  # plot the AUC_score against the ranks of the transcription factors
  
  AUC_plots <- list()
  for (name in names(res)) {
    AUC_plot <- ggplot(data = res[[name]], aes(x = rank, y = AUC_score)) +
      geom_point() +
      geom_line() +
      scale_y_continuous(limits = c(0, 1)) +
      ggtitle(paste("Method:", name))
    AUC_plots[[name]] <- AUC_plot
  }

  # compute the optimal AUC_scores and AUC
  res <- lapply(res, FUN=function(x){
    optimal <- rep(c(TRUE,FALSE),c(length(cofactors),nrow(x)-length(cofactors)))
    cbind(x, optimal)
  })
  maxAUC <- lapply(res, FUN=function(x){
    print(round(sum(cumsum(x$optimal[1:100]) / seq_along(x$optimal[1:100])) / 100,2))
  })
  
  # compute the relative AUC
  relAUC <- list()
  for(method in names(AUC)){
    relAUC[[method]] <- round(AUC[[method]] / maxAUC[[method]], 2)
  }
  
  # return the metrics which were computed above
  return(list(res=res, 
              runtime=runtime, 
              true_rank_table=true_rank_table,
              true_padj=true_padj,
              padj=padj,
              NotTrueTF=NotTrueTF,
              FD=FD,
              AUC=AUC,
              maxAUC=maxAUC,
              relAUC=relAUC,
              rankplot=rankplot,
              padjplot=padjplot,
              AUC_plots=AUC_plots))
}