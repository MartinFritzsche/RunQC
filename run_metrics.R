####### Run Metrics #######

# Master: RunQC_script.R

# Parse runParameters.xml
run_xml <- read_xml(paste0(path, "/runParameters.xml"))


# Generate list to hold high-level run metrics (to be used by report)
run_info <- list(Sequencer = sequencer_type, 
                 RunDate = as.Date(xml_text(xml_find_all(run_xml, ".//RunStartDate")), "%y%m%d"),
                 RunID = xml_text(xml_find_all(run_xml, ".//RunID")),
                 RunNumber = as.integer(xml_text(xml_find_all(run_xml, ".//RunNumber"))),
                 FlowcellSerial = xml_text(xml_find_all(run_xml, ".//SerialNumber"))[1],
                 BufferSerial = xml_text(xml_find_all(run_xml, ".//SerialNumber"))[2],
                 ReagentsSerial = xml_text(xml_find_all(run_xml, ".//SerialNumber"))[3],
                 FlowCellExpiry = as.Date(substr(xml_text(xml_find_all(run_xml, ".//ExpirationDate"))[1], start = 1, stop = 10), "%Y-%m-%d"),
                 BufferExpiry = as.Date(substr(xml_text(xml_find_all(run_xml, ".//ExpirationDate"))[2], start = 1, stop = 10), "%Y-%m-%d"),
                 ReagentsExpiry = as.Date(substr(xml_text(xml_find_all(run_xml, ".//ExpirationDate"))[3], start = 1, stop = 10), "%Y-%m-%d"),
                 Projects = unlist(strsplit(xml_text(xml_find_all(run_xml, ".//ExperimentName")), "_")),
                 Chemistry = xml_text(xml_find_all(run_xml, ".//Chemistry")))

# Generate warning message if any reagent's expiry date is prior to the run date
run_info$ExpiryWarning = run_info$RunDate > run_info$FlowCellExpiry & run_info$RunDate > run_info$BufferExpiry & run_info$RunDate > run_info$ReagentsExpiry

# Fetch info about read design
readDesign_temp <- reads(fc)
readDesign_temp_ls <- list()
for (i in c(1:length(readDesign_temp))) {
  readDesign_temp_ls[[i]] <- data.frame(Number = readDesign_temp[[i]]@number, 
                                        Cycles = readDesign_temp[[i]]@cycles, 
                                        Index = readDesign_temp[[i]]@index)
}
readDesign_temp_df <- do.call("rbind", readDesign_temp_ls)

run_info$RunDesign = ifelse(sum(!readDesign_temp_df$Index) > 1, "Paired-End", "Single-End")
run_info$IndexDesign = ifelse(sum(readDesign_temp_df$Index) > 1, "Dual-Indexed", "Single-Indexed")
run_info$ReadLength = readDesign_temp_df[1, "Cycles"]
run_info$IndexLength = readDesign_temp_df[2, "Cycles"]


# Extract cluster density info
temp <- fc@parsedData$savTileFormat@data
density_temp <- subset(temp, code == 100) # Cluster density
density_temp <- density_temp$value


# Calculate total data generated in Gb
output_temp <- subset(temp, code == 103) # Clusters passing filter
output_temp <- sum(output_temp$value) * sum(!readDesign_temp_df$Index) * run_info$ReadLength / 10^9

run_info$ClusterDensityMean <- round(mean(density_temp / 1000), 0)
run_info$ClusterDensitySD <- round(sd(density_temp / 1000), 0)
run_info$OutputTotal <- output_temp


# Calculate average % >Q30
quality_temp <- fc@parsedData$savQualityFormatV5@data %>%
  select(Q1:Q47) %>%
  colSums() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Quality") %>%
  mutate(Quality = as.integer(substring(Quality, 2))) %>%
  mutate(Threshold = ifelse(Quality > 29, "Pass", "Fail"))

names(quality_temp)[2] <- "Total"

quality_temp <- quality_temp %>%
  mutate(Perc = Total / sum(Total)) %>%
  filter(Threshold == "Pass")

run_info$Q30 <- round(sum(quality_temp$Perc * 100), digits = 2)
run_info$OutputQ30 <- sum(quality_temp$Perc) * run_info$OutputTotal


# Clean-up
rm(run_xml, 
   temp, 
   density_temp, 
   output_temp, 
   readDesign_temp, 
   readDesign_temp_df, 
   readDesign_temp_ls, 
   quality_temp, 
   sequencer_type,
   run_date)
