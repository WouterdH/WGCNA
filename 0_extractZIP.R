#!/usr/bin/env Rscript

### Script that extracts affy .zip files to defined output folder
### Wouter den Hollander 15 Aug 2017

### How to use
##  Rscript 0_extractZIP.R outputDIR zipDIR
##  Where
##    outputDIR must be an existing directory to which the .CEL files are extracted.
##    zipDIR must be a directory from which the affy .zip files are loaded
##  Arguments must be in this order. 


### The actual magic
## Developping / debugging mode
  .dev   <- FALSE
  .debug <- FALSE


## Settings
  options(stringsAsFactors = FALSE, verbose = TRUE)

  if(.dev) {
    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/DM/KIDNEY/', 
    		      '/data/wouter/WGCNA/DATA/anonftp.niehs.nih.gov/drugmatrix/Affymetrix_data/Bulk_data_by_organ/KIDNEY/')
  
    args <- c('/data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/SINGLE/',
              '/data/wouter/WGCNA/DATA/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Kidney/Single/')
  } else {
    args <- commandArgs(trailingOnly = TRUE)
  }
 
  if(length(args) < 2) stop('Need an output directory (outputDIR) and at least one affy zip file input directory (zipDIR_1).\n')

  outputDIR <- args[1]
  inputDIR  <- args[2]


## Pre-run checks
   if(!any(grepl('.zip', list.files(inputDIR)))) stop('No .zip files found in inputDIR.\n')


## Extract .zip files 
  dir.create(outputDIR, recursive = TRUE)
  sys_unzip <- paste0('find ', inputDIR, ' -name "*.zip" -exec unzip {} -d ', outputDIR, ' \\;')
  print(sys_unzip)
  
  system(sys_unzip)


## Take care!
  cat('Thanks for using den Hollander Industries and be sure to contact us again if you need stuff unzipped!\n')

  q()
