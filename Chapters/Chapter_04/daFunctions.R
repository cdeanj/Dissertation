# Author: Chris Dean
# Project: Dissertation Chapter 4
# Date: February 19, 2023
# Purpose: Functions for preprocessing data prior to differential abundance testing
# Website: https://cdeanj.github.io/

library("dplyr")
library("glmmSeq")
library("emmeans")
library("metagenomeSeq")
library("phyloseq")

#' Performs filtering of OTUs from phyloseq object
#'
#' @param psObject An object of class 'phyloseq'
#'
#' @return A filtered phyloseq object
filter_otus <- function(psObj) {
  nOTUs <- nrow(otu_table(psObj))
  nSamples <- ncol(otu_table(psObj))
  
  nOTUsPostFilt <- 0
  
  # Credit to Nicholas Ollberding for filtering threshold
  # https://www.nicholas-ollberding.com/post/introduction-to-phyloseq/
  # Filtering threshold: Retain OTUs with an abundance greater than 10 in at least 10% of samples
  psFilt <- filter_taxa(psObj, function(x) sum(x > 10) > (0.1*length(x)), TRUE)
  psFilt <- prune_taxa(taxa_sums(psFilt) > 0, psFilt)
  
  nOTUsPostFilt <- nrow(otu_table(psFilt))
  
  message("Removed ", nOTUs - nOTUsPostFilt, " OTUs from original ", nOTUs)
  
  return(psFilt)
}

#' Performs filtering of samples from phyloseq object
#'
#' @param psObj An object of class 'phyloseq'
#'
#' @return A filtered phyloseq object
filter_samples <- function(psObj) {
  nSamples <- ncol(otu_table(psObj))
  nSamplesPostFilt <- 0
  
  psFilt <- prune_samples(sample_sums(psObj) > 0, psObj)
  nSamplesPostFilt <- ncol(otu_table(psFilt))
  
  message("Removed ", nSamples - nSamplesPostFilt, " samples from original ", nSamples)
  
  return(psFilt)
}

#' Converts row names in taxonomy table from OTU labels to species labels
#'
#' @param psObj An object of class 'phyloseq'
#'
#' @return A modified phyloseq object
taxa_rownames_to_species <- function(psObj) {
  taxa_names(psObj) <- tax_table(psObj)[,7]
  return(psObj)
}

#' Normalizes count matrix using cumulative sum scaling
#'
#' @param psObj An object of class 'phyloseq'
#'
#' @return A matrix of type double.  Matrix of normalized OTU abundances.
compute_css <- function(psObj) {
  mgs <- phyloseq_to_metagenomeSeq(psObj)
  mgs = cumNorm(mgs, p = 0.5)
  counts <- as(MRcounts(mgs, norm = TRUE, log = FALSE), "matrix")
  return(counts)
}

#' Examines the association between intramammary infection status and the abundance of each OTU
#'
#' @param countMatrix A matrix of type double.  Matrix of OTU abundances where rows = OTUs and cols = samples
#' @param experimentData A table of type list.  Table containing experimental design information about each sample
#'
#' @return An S4 class object of linear mixed effect model results for each OTU
compute_lmms <- function(countData, experimentData) {
  mod <- lmmSeq(~ CaseOrControl + WeeksToInfection + FarmId + Batch + (1 | CowId_Farm), maindata = log(countData + 1), metadata = experimentData, progress = TRUE)
  return(mod)
}

#' Computes estimated marginal means for cases and controls
#'
#' @param mod An S4 class object of linear mixed effect model results for each OTU
#' @param OTUs A vector of OTU names
#'
#' @return A list of two data.frames: (1) contrasts; and (2) emmeans on the response scale
compute_emmeans <- function(mod, otus) {
  contrasts_df <- data.frame()
  emmeans_df <- data.frame()
  
  for(otu in otus) {
    message("Refitting OTU: ", otu)
    fit <- glmmRefit(mod, gene = otu)
    message("Updating reference grid")
    rg <- update(ref_grid(fit), tran = make.tran("genlog", 1))
    message("Computing estimated marginal means")
    emm_obj <- emmeans(regrid(rg), ~ CaseOrControl)
    emm <- emm_obj %>% data.frame()
    message("Computing contrasts")
    contr <- confint(pairs(emm_obj)) %>% data.frame()
    
    contr$otu <- otu
    emm$otu <- otu
    contrasts_df <- rbind(contrasts_df, contr)
    emmeans_df <- rbind(emmeans_df, emm)
  }
  
  return(list(contrasts_df, emmeans_df))
}

#' Computes log2 fold change
#'
#' @param emm An object of type "emmGrid"
#'
#' @return A data frame with log2 fold change between cases and controls
compute_foldchange <- function(emm) {
  cases <- emm %>% filter(CaseOrControl == "Case")
  controls <- emm %>% filter(CaseOrControl == "Control")
  
  nrow <- nrow(cases)
  ncol = 2
  
  ret <- data.frame(matrix(ncol = ncol, nrow = nrow))
  colnames(ret) <- c("otu", "logFC")
  
  if(all(cases$otu == controls$otu)) {
    ret$otu <- cases$otu
    ret$logFC <- log2(cases$emmean) - log2(controls$emmean)
  }
  
  return(ret)
}

#' Returns contrasts between cases and controls on the response scale
#'
#' @param emm A list of data frames containing contrasts and emms
#'
#' @return A data frame of contrasts between cases and controls for each otu
get_contrasts <- function(emm) {
  return(emm[[1]])
}

#' Returns emmeans for cases and controls on the response scale
#'
#' @param emm A list of data frames containing contrasts and emms
#'
#' @return A data frame of emmeans for each otu for cases and controls
get_emmeans <- function(emm) {
  return(emm[[2]])
}