# Author:   Martin Fritzsche
# Date:     21/11/2019
# Purpose:  Develop script-based post-run sequencing QC
# Input:    InterOp folder, runParameters.xml
# Output:   PDF containing plots and tables


####### Commands from savR package #######

# Parse binary files from InterOp folder and create SAV project object
fc <- savR("Test")

# Clusters passing filter per lane
pfBoxplot(fc)

# Get number of sequencing reads (e.g. forward + reverse) (excluding index reads)
directions(fc)

# Get info about reads (F, R, index)
reads(fc)

# Total number of cycles
cycles(fc)

# Info about flow cell (swaths, surfaces, etc)
flowcellLayout(fc)

# Get intensity info as data frame
head(correctedIntensities(fc), n = 10)

# Plot intensities per cycle and tile
plotIntensity(fc)

# Get Quality data per cycle and tile
head(qualityMetrics(fc), n = 10)

# Get encoded tile metrics
head(tileMetrics(fc), n = 4)

# 
head(extractionMetrics(fc), n = 2)

# Extract Run ID
run(fc)


####### Set Up #######

source(file = "libraries.R")

# [ TO DO ] Use that as input for script
sequencer_type <- "NextSeq"
run_date <- "191106"

path <- file.path(paste("./data", sequencer_type, run_date, sep = "/"))

fc <- savR(path) # Generate SAV project object from Illumina binaries

#loadfonts(device = "win")

####### Run Metrics #######

run_xml <- read_xml(paste0(path, "/runParameters.xml"))

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

run_info$ExpiryWarning = run_info$RunDate > run_info$FlowCellExpiry & run_info$RunDate > run_info$BufferExpiry & run_info$RunDate > run_info$ReagentsExpiry

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

temp <- fc@parsedData$savTileFormat@data
density_temp <- subset(temp, code == 100) # Cluster density
density_temp <- density_temp$value

output_temp <- subset(temp, code == 103) # Clusters passing filter
output_temp <- sum(output_temp$value) * sum(!readDesign_temp_df$Index) * run_info$ReadLength / 10^9

run_info$ClusterDensityMean <- round(mean(density_temp / 1000), 0)
run_info$ClusterDensitySD <- round(sd(density_temp / 1000), 0)
run_info$OutputTotal <- output_temp

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


rm(run_xml, temp, density_temp, output_temp, readDesign_temp, readDesign_temp_df, readDesign_temp_ls, quality_temp)
  
####### Quality Metric ######

if (run_info$Sequencer == "MiSeq") {
  quality_df <- qualityMetrics(fc) # This does not work for NextSeq!
} else {
  quality_df <- fc@parsedData$savQualityFormatV5@data
}



#QScore Distribution Overview
quality_overview <- quality_df %>%
  select(Q1:Q47)

quality_overview <- as.data.frame(colSums(quality_overview))
quality_overview <- tibble::rownames_to_column(quality_overview, "Quality")
quality_overview <- quality_overview %>%
  mutate(Quality = as.integer(substring(Quality, 2))) %>%
  mutate(Threshold = ifelse(Quality > 29, "Pass", "Fail"))
  
names(quality_overview)[2] <- "Total"

quality_overview <- quality_overview %>%
  mutate(Perc = round(Total / sum(Total) * 100, digits = 2)) 

quality_perc <- quality_overview %>%
  group_by(Threshold) %>%
  summarise(Sum = sum(Total)) %>%
  mutate(Percent = Sum / sum(Sum))
quality_perc <- as.double(quality_perc[2, "Percent"])

ggplot(quality_overview, aes(x = as.factor(Quality), y = Total / 10^9, fill = Threshold)) +
  geom_bar(stat = "identity") +
  annotate("text", label = paste(">= Q30", scales::percent(quality_perc, accuracy = 0.01),  sep = "\n"), x = -Inf, y = Inf, size = 8, hjust = -6, vjust = 2.5, colour = "chartreuse4") +
  geom_text(aes(label = Perc), size = 2.8, vjust = -0.5, colour = "grey50") +
  xlab("Q Score") +
  ylab("Total Bases (Billions)") +
  ggtitle("Q Score Distribution") +
  theme_minimal() +
  theme(legend.position = "none")


# Quality overview by cycle
quality_overview_cycle <- quality_df %>%
  pivot_longer(-c(lane, cycle, tile), names_to = "Quality", values_to = "Value") %>%
  mutate(Quality = as.integer(substring(Quality, 2))) %>%
  mutate(Threshold = ifelse(Quality > 29, "Pass", "Fail")) %>%
  group_by(lane, cycle, tile, Threshold) %>%
  summarise(Sum = sum(Value)) %>%
  pivot_wider(names_from = "Threshold", values_from = "Sum") %>%
  mutate(Perc_Pass = Pass / (Pass + Fail) * 100)

read_intercepts <- c(R1 = 1, 
                     R2 = 1 + run_info$ReadLength, 
                     R3 = 1 + run_info$ReadLength + run_info$IndexLength, 
                     R4 = 1 + run_info$ReadLength + 2 * run_info$IndexLength)

ggplot(quality_overview_cycle, aes(x = cycle, y = Perc_Pass, group = cycle)) + 
  geom_boxplot(outlier.size = 1, outlier.colour = "red") +
  geom_vline(xintercept = read_intercepts, colour = "blue", linetype = "dashed") +
  annotate("text", x = read_intercepts, y = 104, label = names(read_intercepts), colour = "blue", size = 3.6, hjust = -0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face = "bold")) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 10), expand = c(0, 0)) +
  ylab("% Clusters >= Q30") +
  ggtitle("Quality Overview by Cycle")



# Quality overview by tile and cycle
quality_overview_tile <- select(quality_overview_cycle, -c(Fail, Pass))

r1  <- subset(quality_overview_tile, cycle <= run_info$ReadLength)
r23 <- subset(quality_overview_tile, cycle > run_info$ReadLength & cycle <= run_info$ReadLength + 2 * run_info$IndexLength)
r4  <- subset(quality_overview_tile, cycle > run_info$ReadLength + 2 * run_info$IndexLength)

draw_tile_plot <- function(df) {
  ggplot(df, aes(x = cycle, y = as.factor(tile), fill = Perc_Pass)) +
    geom_tile(colour = "grey33") +
    theme_minimal() +
    theme(legend.position = "top",
          legend.justification = 0,
          legend.title = element_text(face = "bold"),
          axis.ticks = element_line(),
          axis.text.y = element_text(size = 4),
          plot.title = element_text(face = "bold")) +
    scale_fill_gradient(low = "red", high = "blue", name = "% >= Q30") +
    scale_y_discrete("Tile", expand = c(0, 0))
} 

draw_tile_plot(subset(r1, lane == 1)) + 
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = run_info$ReadLength, by = 10), expand = c(0, 0)) +
  ggtitle("R1 Tile Quality - Lane 1")

draw_tile_plot(subset(r1, lane == 2)) + 
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = run_info$ReadLength, by = 10), expand = c(0, 0)) +
  ggtitle("R1 Tile Quality - Lane 2")

draw_tile_plot(subset(r1, lane == 3)) + 
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = run_info$ReadLength, by = 10), expand = c(0, 0)) +
  ggtitle("R1 Tile Quality - Lane 3")

draw_tile_plot(subset(r1, lane == 4)) + 
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = run_info$ReadLength, by = 10), expand = c(0, 0)) +
  ggtitle("R1 Tile Quality - Lane 4")


draw_tile_plot(subset(r23, lane == 1)) + 
  scale_x_continuous("Cycle", breaks = seq(from = run_info$ReadLength + 1, to = run_info$ReadLength + 1 + 2 * run_info$IndexLength, by = 1), expand = c(0, 0)) +
  ggtitle("Index Reads Tile Quality - Lane 1")

draw_tile_plot(subset(r23, lane == 2)) + 
  scale_x_continuous("Cycle", breaks = seq(from = run_info$ReadLength + 1, to = run_info$ReadLength + 1 + 2 * run_info$IndexLength, by = 1), expand = c(0, 0)) +
  ggtitle("Index Reads Tile Quality - Lane 2")

draw_tile_plot(subset(r23, lane == 3)) + 
  scale_x_continuous("Cycle", breaks = seq(from = run_info$ReadLength + 1, to = run_info$ReadLength + 1 + 2 * run_info$IndexLength, by = 1), expand = c(0, 0)) +
  ggtitle("Index Reads Tile Quality - Lane 3")

draw_tile_plot(subset(r23, lane == 4)) + 
  scale_x_continuous("Cycle", breaks = seq(from = run_info$ReadLength + 1, to = run_info$ReadLength + 1 + 2 * run_info$IndexLength, by = 1), expand = c(0, 0)) +
  ggtitle("Index Reads Tile Quality - Lane 4")


draw_tile_plot(subset(r4, lane == 1)) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 2 * run_info$ReadLength + 2 * run_info$IndexLength, by = 10), expand = c(0, 0)) +
  ggtitle("R4 Tile Quality - Lane 1")

draw_tile_plot(subset(r4, lane == 2)) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 2 * run_info$ReadLength + 2 * run_info$IndexLength, by = 10), expand = c(0, 0)) +
  ggtitle("R4 Tile Quality - Lane 2")

draw_tile_plot(subset(r4, lane == 3)) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 2 * run_info$ReadLength + 2 * run_info$IndexLength, by = 10), expand = c(0, 0)) +
  ggtitle("R4 Tile Quality - Lane 3")

draw_tile_plot(subset(r4, lane == 4)) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 2 * run_info$ReadLength + 2 * run_info$IndexLength, by = 10), expand = c(0, 0)) +
  ggtitle("R4 Tile Quality - Lane 4")


# Quality by cycle (separate bins)
quality_cycle <- quality_df %>%
  group_by(cycle) %>%
  summarise_at(vars(Q1:Q47), mean) %>%
  pivot_longer(-c(cycle), names_to = "Quality", values_to = "Value") %>%
  mutate(quality = as.integer(substring(Quality, 2)))

ggplot(quality_cycle, aes(x = cycle, y = Value)) +
  geom_line() + 
  facet_wrap(vars(as.factor(Quality)))



# Quality by cycle (separate bins)
# quality_tile <- quality_df %>%
#   pivot_longer(-c(lane, tile, cycle), names_to = "Quality", values_to = "Value") %>%
#   group_by(tile, Quality) %>%
#   summarise(Mean = mean(Value), SD = sd(Value)) %>%
#   mutate(Surface = ifelse(tile < 2000, "Top", "Bottom")) %>%
#   mutate(Quality = as.integer(substring(Quality, 2)))
# 
# ggplot(quality_tile, aes(x = as.factor(tile), y = Mean, fill = Surface)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.1) +
#   facet_wrap(vars(as.factor(Quality))) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("Tile") +
#   ylab("Clusters")

# quality_tile <- quality_df %>%
#   pivot_longer(-c(lane, tile, cycle), names_to = "quality", values_to = "value") %>%
#   mutate(surface = ifelse(tile < 2000, "Top", "Bottom")) %>%
#   mutate(quality = as.integer(substring(quality, 2)))
# 
# ggplot(quality_tile, aes(x = as.factor(tile), y = value, colour = surface )) +
#   geom_point(alpha = 0.5) +
#   facet_wrap(vars(as.factor(quality))) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("Tile") +
#   ylab("Clusters") +
#   ggtitle("Cluster Quality per Tile")



# Quality per tile and cycle
quality_tile_cycle <- quality_df %>%
  pivot_longer(-c(tile, cycle, lane), names_to = "Quality", values_to = "Value") %>%
  mutate(Quality = as.integer(substring(Quality, 2))) %>%
  mutate(Tile = as.integer(tile)) %>%
  filter(Quality %in% c(14, 21, 27, 32, 36)) 
  
  
ggplot(quality_tile_cycle, aes(x = cycle, y = Value / 1000, colour = Tile)) +
  geom_line() + 
  facet_wrap(vars(as.factor(Quality))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 20)) +
  ylab("k Clusters") +
  ggtitle("Cluster Quality per Cycle")


# Cleanup
rm(list = c("draw_tile_plot", "quality_df", "quality_overview", "quality_overview_cycle", "quality_overview_tile", "quality_perc", "quality_cycle", "quality_tile_cycle", "r1", "r23", "r4"))



####### Tile Metrics #######

tiles_df <- tileMetrics(fc)

# Translate codes into variable names
# NOTE: Phasing and prephasing codes depend on number of sequencing reads
tiles_df$code <- recode(tiles_df$code, 
                        "100" = "Density", 
                        "101" = "Density_PF", 
                        "102" = "Clusters", 
                        "103" = "Clusters_PF", 
                        "200" = "Phasing_R1", 
                        "201" = "Prephasing_R1",
                        "202" = "Phasing_R2",
                        "203" = "Prephasing_R2",
                        "204" = "Phasing_R3",
                        "205" = "Prephasing_R4",
                        "206" = "Phasing_R4",
                        "207" = "Prephasing_R4",
                        "300" = "PercAligned_R1", 
                        "303" = "PercAligned_R4",
                        "400" = "Control_Lane")

# Keep table in long format, separate tiles into top and bottom surface
tiles_df_long <- tiles_df %>%
  mutate(tile = as.character(tile)) %>%
  mutate(Surface = ifelse(gsub("^(.)(.)(.)(.)(.)$", "\\1", tile) == "1", "Top", "Bottom")) %>%
  mutate(Swath = gsub("^(.)(.)(.)(.)(.)$", "\\2", tile)) %>%
  mutate(Camera = gsub("^(.)(.)(.)(.)(.)$", "\\3", tile)) %>%
  mutate(tile = substring(tile, 4))


# Cluster Plot
ggplot(subset(tiles_df_long, code == "Clusters_PF"), aes(x = Swath, y = as.factor(tile), fill = value / 10^6)) +
  geom_tile() +
  facet_grid(rows = vars(Surface), cols = vars(Camera)) +
  theme_minimal() +
  ggtitle("Clusters per Camera") +
  ylab("Tile") +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))






# Phasing Plot
# ggplot(subset(tiles_df_long, code %in% c("Phasing_R1", "Phasing_R4", "Prephasing_R1", "Prephasing_R4")), aes(x = Surface, y = as.factor(tile), fill = value * 100)) +
#   geom_tile() + 
#   facet_wrap(vars(code)) +
#   theme_minimal() +
#   scale_fill_distiller(palette = "Set3", name = "Value * 10^6") +
#   ggtitle("Flow Cell Info - Phasing") +
#   ylab("Tile") +
#   theme(strip.text.x = element_text(size = 12, face = "bold"),
#         axis.title.x = element_text(size = 14),
#         axis.title.y = element_text(size = 14))

# Cleanup
rm(list = c("tiles_df", "tiles_df_long"))

####### Intensity Metrics #######

intensity_df <- correctedIntensities(fc)

intensity_df <- intensity_df %>%
  mutate(Surface = ifelse(tile < 2000, "Top", "Bottom")) %>%
  mutate(Per_A = num_A / (num_A + num_C + num_G + num_T) * 100) %>%
  mutate(Per_C = num_C / (num_A + num_C + num_G + num_T) * 100) %>%
  mutate(Per_G = num_G / (num_A + num_C + num_G + num_T) * 100) %>%
  mutate(Per_T = num_T / (num_A + num_C + num_G + num_T) * 100)

base_tile_avg <- intensity_df %>%
  select(tile, cycle, Per_A, Per_C, Per_G, Per_T, Surface) %>%
  group_by(cycle, Surface) %>%
  summarise(A = mean(Per_A), 
            C = mean(Per_C), 
            G = mean(Per_G),
            `T` = mean(Per_T)) %>%
  pivot_longer(-c(cycle, Surface), names_to = "Base", values_to = "Value")

ggplot(base_tile_avg, aes(x = cycle, y = Value, colour = Base)) +
  geom_line() +
  geom_vline(xintercept = read_intercepts, colour = "blue", linetype = "dashed") +
  annotate("text", x = read_intercepts, y = 104, label = names(read_intercepts), colour = "blue", size = 3.6, hjust = -0.3) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold")) +
  ylab("% Base") +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 10), expand = c(0, 0)) +
  ggtitle("Base Composition per Cycle")



base_tile_sd <- intensity_df %>%
  select(tile, cycle, Per_A, Per_C, Per_G, Per_T) %>%
  group_by(cycle) %>%
  summarise(A = sd(Per_A), 
            C = sd(Per_C), 
            G = sd(Per_G),
            `T` = sd(Per_T)) %>%
  pivot_longer(-c(cycle), names_to = "Base", values_to = "SD")

ggplot(base_tile_sd, aes(x = cycle, y = SD, colour = Base)) +
  geom_line() +
  theme_minimal() +
  geom_vline(xintercept = read_intercepts, colour = "blue", linetype = "dashed") +
  annotate("text", x = read_intercepts, y = 2, label = names(read_intercepts), colour = "blue", size = 3.6, hjust = -0.3) +
  ylab("% Base Standard Deviation") +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 10), expand = c(0, 0)) +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold")) +
  ggtitle("Base Call Dispersion across Tiles")



base_tile_sep <- intensity_df %>%
  select(tile, cycle, Per_A, Per_C, Per_G, Per_T) %>%
  pivot_longer(-c(cycle, tile), names_to = "Base", values_to = "Percent")

ggplot(base_tile_sep, aes(x = cycle, y = as.factor(tile), fill = Percent)) +
  geom_raster() +
  scale_fill_distiller(palette = "Blues") +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 10), expand = c(0, 0)) +
  facet_wrap(vars(Base)) +
  theme_minimal() +
  ggtitle("Tile Base Calls") +
  ylab("Tile") +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))



# Calculate Shannon entropy for each position (log basis is 2)
shannon_entropy <- function(input) {
  results <- numeric(length = 4)
  i <- 1
  for (value in input) {
    value <- -value * log2(value)
    if (is.nan(value)) value <- 0
    results[i] <- value
    i <- i + 1 
  }
  results <- sum(results)
  return(results)
}

base_tile_entropy <- intensity_df %>%
  select(tile, cycle, Per_A, Per_C, Per_G, Per_T) %>%
  mutate_at(c("Per_A", "Per_C", "Per_G", "Per_T"), function(x) x / 100)

base_tile_temp <- select(base_tile_entropy, -c(tile, cycle))
base_tile_temp <- mutate(base_tile_temp, Shannon = apply(base_tile_temp, 1, shannon_entropy))
base_tile_entropy$Shannon <- base_tile_temp$Shannon
rm(base_tile_temp)

# ggplot(base_tile_entropy, aes(x = cycle, y = as.factor(tile), fill = Shannon)) +
#   geom_raster() +
#   scale_fill_distiller(palette = "Blues") +
#   scale_x_continuous("Cycle", breaks = seq(from = 0, to = 520, by = 10), expand = c(0, 0)) +
#   theme_minimal() +
#   ggtitle("Base Entropy per Cycle") +
#   ylab("Tile") +
#   theme(legend.position = "top",
#         legend.justification = 0,
#         legend.title = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(subset(base_tile_entropy, cycle <= 251), aes(x = cycle, y = as.factor(tile), fill = Shannon)) +
  geom_tile(colour = "grey33") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold"),
        axis.ticks = element_line()) +
  scale_fill_gradient(low = "red", high = "blue", name = "Shannon Entropy") +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = run_info$ReadLength, by = 10), expand = c(0, 0)) +
  scale_y_discrete("Tile", expand = c(0, 0)) + 
  ggtitle("Forward Read - Base Entropy per Cycle")

ggplot(subset(base_tile_entropy, cycle > 251 & cycle <= 267), aes(x = cycle, y = as.factor(tile), fill = Shannon)) +
  geom_tile(colour = "grey33") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold"),
        axis.ticks = element_line()) +
  scale_fill_gradient(low = "red", high = "blue", name = "Shannon Entropy") +
  scale_x_continuous("Cycle", breaks = seq(from = run_info$ReadLength + 1, to = run_info$ReadLength + 1 + 2 * run_info$IndexLength, by = 1), expand = c(0, 0)) +
  scale_y_discrete("Tile", expand = c(0, 0)) + 
  ggtitle("Index Reads - Base Entropy per Cycle")

ggplot(subset(base_tile_entropy, cycle > 267), aes(x = cycle, y = as.factor(tile), fill = Shannon)) +
  geom_tile(colour = "grey33") +
  theme_minimal() +
  theme(legend.position = "top",
        legend.justification = 0,
        legend.title = element_text(face = "bold"),
        axis.ticks = element_line()) +
  scale_fill_gradient(low = "red", high = "blue", name = "Shannon Entropy") +
  scale_x_continuous("Cycle", breaks = seq(from = 0, to = 2 * run_info$ReadLength + 2 * run_info$IndexLength, by = 10), expand = c(0, 0)) +
  scale_y_discrete("Tile", expand = c(0, 0)) + 
  ggtitle("Reverse Read - Base Entropy per Cycle")
