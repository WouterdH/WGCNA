#!/usr/bin/env Rscript

### Script that loads supplied rmaNormalized .RData file, then calculates the FC and log2 transform those
### for each compound * timepoint * concentration combination vs respecitve control treatment.

### The script is currently tailored towards seperate use of DM and TGG TXG data.

### Wouter den Hollander 29 Sept 2017

### How to use
##  Rscript 2_selectProbes_calculateFC_scaleFC.R outputDIR rmaFile.RData 
##  Where
##    outputDIR         must be a directory to which the scaled FC  data is saved.
##    rmaFile.RData     RMA normalized TXG data from either TGG or DM.
##                      Since the script is tailored towards DM & TGG, contrasts will
##                      be clear from this file.
##  Arguments must be in this order. 


### The actual magic
## Developping / debugging mode
  .dev   <- FALSE
  .debug <- FALSE


## Functions
  parseTG <- function(rmaObject = rmaData) {
    # This functions adds a column to pData(rmaData) that reflects a unique combination
    # of compound, time, dose and single/repeat dosing for TG data.

    EXP_UID <- toupper(gsub(' ', '', apply(pData(rmaObject)[, c('COMPOUND ABBR.', 'SACRI_PERIOD', 'DOSE_LEVEL', 'SIN_REP_TYPE')], 
                                       1, 
                                       paste, 
                                       collapse = '_')))

    pData(rmaObject)$EXP_UID <- EXP_UID
    return(rmaObject)
  }


  fitLimmaSubset <- function(EXP_UID, rmaObject = rmaData, returnFit = FALSE, verbose = TRUE) {
    if(.debug) {
      EXP_UID <- 'APAP_4DAY_LOW_REPEAT'
      rmaObject    <- rmaData
    }

    if(verbose) cat(EXP_UID, '\n')

    EXP_UID_ctrl <- grep('CONTROL', sapply(c('LOW', 'MIDDLE', 'HIGH'), function(.dose) gsub(.dose, 'CONTROL', EXP_UID) ), value = TRUE)

    .rmaObject <- rmaObject[, which(pData(rmaObject)$EXP_UID == EXP_UID | pData(rmaObject)$EXP_UID == EXP_UID_ctrl)]
      pData(.rmaObject)$CONTRAST <- as.numeric(pData(.rmaObject)$DOSE_FACTOR != 'CTRL')

    design <- model.matrix(~ 1 + CONTRAST, data = pData(.rmaObject))
    fit <- eBayes(lmFit(.rmaObject, design))

    if(returnFit) return(fit) else return( data.frame(topTable(fit, coef = 2, number = nrow(.rmaObject)), experiment = EXP_UID) )
  }

## Settings
  options(stringsAsFactors = FALSE, verbose = TRUE)

  if(.dev) {
    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/KIDNEY/RMA/',
              '/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/KIDNEY/RMA/RMA_2017-09-28.RData')

    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/RMA/',
              '/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/RMA/RMA_2017-09-28.RData')
  } else {
    args <- commandArgs(trailingOnly = TRUE)
  }

  if(length(args) < 2) stop('Need outputdir and rmaFile to work with.\n')

  outputDIR <- args[1]
  rmaFile   <- args[2]

  if(!file.exists(rmaFile)) stop('Supplied rmaFile does not seem to exist.\n')


## Required packages
  library(affy)
  library(panp)
  library(parallel)
  library(limma)


## Load & parse RMA data
  load(rmaFile)
  if(!exists('rmaData')) stop('Supplied rmaFile does not contain rmaData object.\n')

  if(all(pData(rmaData)$SOURCE == 'TG')) rmaData <- parseTG()
  if(all(pData(rmaData)$SOURCE == 'DM')) rmaData <- parseDM()

  EXP_UIDs <- unique(grep('CONTROL', pData(rmaData)$EXP_UID, value = TRUE, invert = TRUE))

## DEG analysis
  listDEG <- mclapply(EXP_UIDs, function(.EXP_UID) fitLimmaSubset(EXP_UID = .EXP_UID), mc.cores = 24)
  listDEG <- mclapply(listDEG, function(.deg) .deg[order(rownames(.deg)), ], mc.cores = 24 )
    names(listDEG) <- sapply(listDEG, function(.deg) .deg$experiment[1])

  dfDEG <- data.frame(lapply(listDEG, function(.deg) .deg[, 'logFC', drop = FALSE]))
    colnames(dfDEG) <- names(listDEG)


## Z-scaling of genes
  dfDEG_z <- apply(dfDEG, 1, function(x) (x - mean(x)) / sd(x) )


## Save data
  dir.create(outputDIR, recursive = TRUE)
  save(dfDEG_z, file = paste0(outputDIR, 'zScaledFCs_RMA_', Sys.Date(), '.RData'))


## Take care!
  cat('Thanks for using den Hollander Industries and be sure to contact us again if you need FC values calculated & scaled!\n')

  q()














