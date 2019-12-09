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

source(file = "run_metrics.R")

####### Quality Metrics #######

source(file = quality_metrics.R)

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
