# Author: Chris Dean
# Project: Dissertation Chapter 4
# Date: February 19, 2023
# Purpose: Differential abundance testing
# Website: https://cdeanj.github.io/

library("dplyr")
library("ggplot2")

source("daFunctions.R")

path.rds <- "../RDS"

ps <- readRDS(file.path(path.rds, "ps_shotgun_02192023.rds"))
ps

psPrune <- filter_otus(ps) # filter OTUs based on abundance and prevalence
psPrune <- filter_samples(psPrune) # remove samples with no OTUs
psPrune <- taxa_rownames_to_species(psPrune)

counts <- compute_css(psPrune) # normalize counts using css
metadata <- sample_data(psPrune) %>% data.frame() # extract meta data

mod <- compute_lmms(counts, metadata) # build lmms for each otu
mod <- glmmQvals(mod) # calculate q values

modSummary <- summary(mod) %>% data.frame() # generate model summaries
features <- row.names(modSummary) # extract otu names

emmList <- compute_emmeans(mod, features) # compute estimated marginal means for each lmm in mod

contrastsDF <- get_contrasts(emmList) # return contrasts between cases and controls
emmeansDF <- get_emmeans(emmList) # return emmeans for cases and controls

foldchangeDF <- compute_foldchange(emmeansDF) # compute log2 fold change

saveRDS(psPrune, file.path(path.rds, "da_ps_object.rds"))
saveRDS(modSummary, file.path(path.rds, "da_mod_object.rds"))
saveRDS(contrastsDF, file.path(path.rds, "da_contrasts_object.rds"))
saveRDS(emmeansDF, file.path(path.rds, "da_emm_object.rds"))
saveRDS(foldchangeDF, file.path(path.rds, "da_log2fc_object.rds"))
