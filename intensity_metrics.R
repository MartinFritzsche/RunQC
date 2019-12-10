####### Intensity Metrics #######

# Master: RunQC_script.R


# Parse intensity metrics
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