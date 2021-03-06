#!/usr/bin/env Rscript

### Script that reads .CEL files and generates an RMA normalized .RData file in the 
### defined outputDIR

### Wouter den Hollander 15 Aug 2017

### How to use
##  Rscript 1_readCEL_runRMA.R outputDIR customCDF celDIR_1 celDIR_2 ... celDIR_n
##  Where
##    outputDIR   must be a directory to which the normalized expression data is saved.
##    customCDF   must be an R PACKAGE INSTALLED customCDF file
##    celDIR_n    must be a directory that contains subdirectories with .CEL files and annotation (DM)	  
##  Arguments must be in this order. 


### The actual magic
## Developping / debugging mode
  .dev   <- FALSE
  .debug <- FALSE


## Functions
  .mgsub <- function(pattern, replacement, x, ...) {
    n <- length(pattern)
    if (n != length(replacement)) stop("pattern and replacement do not have the same length.")
    
    rslt <- x
    for (i in 1:n) rslt[grep(pattern[i], x, ...)] = replacement[i]
    
    return(rslt)
  }

  readAnnotation <- function(ann_file) {
    ## Function is tailored towards IN_VIVO rat data from TG & DM.

    if(.debug) {
      ann_file <- '/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/SINGLE//indomethacin.Rat.in_vivo.Kidney.Single/Attribute.tsv' # TG
      ann_file <- '/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/Affymetrix_KIDNEY_A/ATORVASTATIN/ATORVASTATIN'               # DM
    }

    cat(ann_file, '\n')

    ann_data <- data.frame(do.call(rbind, lapply(readLines(ann_file), function(line) unlist(strsplit(gsub('\t', '___', line), '___') )))) 
    ann_data <- ann_data[ann_data[,2] != 'NA', ]
      colnames(ann_data) <- toupper(ann_data[1, ])
      rownames(ann_data) <- paste0(ann_data[, 1], '.CEL')
    ann_data <- ann_data[2:nrow(ann_data), ]  

    if(!any(colnames(ann_data) == 'DOSE_UNIT')) {   # Apparently a DM annotation file
      .various_dose_ind <- which(ann_data$DOSE == 'various')
      .defined_dose_ind <- which(ann_data$DOSE != 'various')

      .dose_data <- do.call(rbind, strsplit(ann_data$DOSE[.defined_dose_ind], ' '))[, 1]
      .dose_unit <- unique(do.call(rbind, strsplit(ann_data$DOSE[.defined_dose_ind], ' '))[, 2])

      if(length(.dose_unit) != 1) stop('Dose unit is not of a single type.\n')

      ann_data$DOSE[.various_dose_ind] <- 0
      ann_data$DOSE[.defined_dose_ind] <- .dose_data
      ann_data$DOSE_UNIT <- .dose_unit

      colnames(ann_data)[which(colnames(ann_data) == 'CHEMICAL')] <- 'COMPOUND_NAME'
      colnames(ann_data)[which(colnames(ann_data) == 'ORGAN')] <- 'ORGAN_ID'
      colnames(ann_data)[which(colnames(ann_data) == 'ARRAY')] <- 'BARCODE'  
        ann_data$BARCODE <- gsub('.CEL', '', ann_data$BARCODE) 

      if(any(ann_data$DOSE == 0)) {
        ann_data$DOSE_FACTOR <- c('CTRL', 'LOW', 'MED', 'HI')[as.numeric(factor(as.numeric(ann_data$DOSE), ordered = TRUE))]
      } else {
        ann_data$DOSE_FACTOR <- c('LOW', 'MED', 'HI')[as.numeric(factor(as.numeric(ann_data$DOSE), ordered = TRUE))]
      }

      ann_data$DATA_SOURCE <- 'DM'
    } else {  # Apparently a TG annotation file
      ann_data$DOSE_FACTOR <- .mgsub(c('CONTROL', 'LOW', 'MIDDLE', 'HIGH'), c('CTRL', 'LOW', 'MED', 'HI'), toupper(ann_data$DOSE_LEVEL))
      ann_data$DATA_SOURCE <- 'TG'
      ann_data$TIME <- factor(as.numeric(gsub(' day', '', gsub(' hr', '', ann_data$SACRI_PERIOD)), ordered = TRUE))
    }

    return(ann_data)
  }


## Settings
  options(stringsAsFactors = FALSE, verbose = TRUE)

  if(.dev) {
    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/KIDNEY/RMA/',
              'rat2302rnentrezgcdf',
              '/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/KIDNEY/')

    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/RMA/',
              'rat2302rnentrezgcdf',
              '/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/')
  } else {
    args <- commandArgs(trailingOnly = TRUE)
  }
 
  if(length(args) < 3) stop('Need outputdir, customCDF and at least one celDIR.\n')

  outputDIR <- args[1]
  customCDF <- args[2]
  celDIR    <- args[3:length(args)]
  

## Required packages
  library(affy)
  library(rat2302rnentrezgcdf)
  library(plyr)


## Define files
  celDIR_files <- list.files(celDIR, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    if(!any(grepl('.CEL', celDIR_files))) stop('Could not find .CEL files.\n')
  
  cel_files <- celDIR_files[ grepl('.CEL', celDIR_files)]
  ann_files <- celDIR_files[!grepl('.CEL', celDIR_files)]

  cat(paste0('Found ', length(cel_files), ' CEL files to RMA.\n'))
  cat(paste0('Found ', length(ann_files), ' annotation files along with them.\n'))

  cel_files_unique <- unique(sapply(cel_files, function(x) grep('CEL', unlist(strsplit(x, '/')), value = T )))
  if(length(cel_files) != length(cel_files_unique)) { # This might fuck-up TG data; not tested.
    cat('Non-unique .CEL files found. Dropping duplicates.\n')

    cel_files <- cel_files[which(!duplicated(sapply(cel_files, function(x) grep('CEL', unlist(strsplit(x, '/')), value = T ))))]
  }

## RMA
  rmaData  <- justRMA(filenames = cel_files, celfile.path = '', cdfname = customCDF)
    sampleData <- do.call(rbind.fill, lapply(ann_files, readAnnotation))
    sampleData <- sampleData[which(!duplicated(sampleData$BARCODE)), ]
      rownames(sampleData) <- paste0(sampleData$BARCODE, '.CEL')
    pData(rmaData) <- sampleData
    exprs(rmaData) <- exprs(rmaData)[, rownames(pData(rmaData))]


## Save data
  dir.create(outputDIR, recursive = TRUE)
  save(rmaData, file = paste0(outputDIR, 'RMA_', Sys.Date(), '.RData'))


## Take care!
  cat('Thanks for using den Hollander Industries and be sure to contact us again if you need stuff RMA normalized!\n')

  q()





  