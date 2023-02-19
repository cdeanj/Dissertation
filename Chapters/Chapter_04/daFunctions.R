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

#' Computes estimated marginal means for each OTU between cases and controls
#'
#' @param mod An S4 class object of linear mixed effect model results for each OTU
#' @param OTUs A vector of OTU names
#'
#' @return A data frame of marginal contrasts between cases and controls for each OTU
compute_contrasts <- function(mod, otus) {
  df <- data.frame()
  
  for(otu in otus) {
    message("Refitting OTU: ", otu)
    fit <- glmmRefit(mod, gene = otu)
    message("Updating reference grid")
    rg <- update(ref_grid(fit), tran = make.tran("genlog", 1))
    message("Computing estimated marginal means")
    emm <- emmeans(regrid(rg), ~ CaseOrControl)
    message("Computing contrasts")
    contr <- confint(pairs(emm)) %>% data.frame()
    
    contr$otu <- otu
    df <- rbind(df, contr)
  }
  
  return(df)
}


