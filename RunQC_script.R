# Author:   Martin Fritzsche
# Date:     21/11/2019
# Purpose:  Develop script-based post-run sequencing QC
# Input:    InterOp folder, 
#           runParameters.xml
#           SampleSheet.csv
# Output:   HTML rendered markdown report containing plots and tables



####### Set Up #######

source(file = "libraries.R")

# [ TO DO ] Use that as input for script
sequencer_type <- "NextSeq"
run_date <- "200122"

path <- file.path(paste("./data", sequencer_type, run_date, sep = "/"))

fc <- savR(path) # Generate SAV project object from Illumina binaries

#loadfonts(device = "win")

####### Run Metrics #######

source(file = "run_metrics.R")

####### Quality Metrics #######

source(file = "quality_metrics.R")

####### Tile Metrics #######

source(file = "tile_metrics.R")

####### Intensity Metrics #######

source(file = "intensity_metrics.R")
