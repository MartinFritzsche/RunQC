####### Tile Metrics #######

# Master: RunQC_script.R

# Parse tile metrics
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
