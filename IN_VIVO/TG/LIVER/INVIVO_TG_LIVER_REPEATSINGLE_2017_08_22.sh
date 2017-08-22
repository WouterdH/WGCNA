# Rscript ~/gitProjects/WGCNA/0_extractZIP.R /data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/SINGLE/ /data/wouter/WGCNA/DATA/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Single/
# Rscript ~/gitProjects/WGCNA/0_extractZIP.R /data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/REPEAT/ /data/wouter/WGCNA/DATA/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Repeat/

Rscript ~/gitProjects/WGCNA/1_readCEL_runRMA.R /data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/ rat2302rnentrezgcdf /data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/SINGLE/ /data/wouter/WGCNA/OUTPUT/IN_VIVO/TG/LIVER/REPEAT/ 

