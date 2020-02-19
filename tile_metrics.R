####### Tile Metrics #######

# Master: RunQC_script.R
# Requires Upstream: RunQC_script.R


###### Parse tile metrics ######

tiles_df <- tileMetrics(fc)

###### Translate codes into variable names ######

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

if (run_info$Sequencer == "NextSeq") {
  tiles_df_long <- tiles_df %>%
    mutate(Surface = ifelse(substring(as.character(tile), first = 1, last = 1) == "1", "Top", "Bottom")) %>%
    mutate(Swath = substring(as.character(tile), first = 2, last = 2)) %>%
    mutate(Camera = substring(as.character(tile), first = 3, last = 3))
  tiles_df_long$Swath <- factor(tiles_df_long$Swath, levels = c("1", "2", "3"), labels = c("Swath 1", "Swath 2", "Swath 3"))
  tiles_df_long$Camera <- factor(tiles_df_long$Camera, levels = c("1", "2", "3", "4", "5", "6"), labels = c("Camera 1", "Camera 2", "Camera 3", "Camera 4", "Camera 5", "Camera 6"))
} else {
  tiles_df_long <- tiles_df %>%
    mutate(Surface = ifelse(substring(as.character(tile), first = 1, last = 1) == "1", "Top", "Bottom")) %>%
    mutate(tile = substring(tile, 3))
} 



# Cluster Plot MiSeq
if (run_info$Sequencer == "MiSeq") {
  ggplot(subset(tiles_df_long, code %in% c("Density", "Density_PF")), aes(x = as.factor(lane), y = as.factor(tile), fill = value)) +
    geom_tile() +
    facet_grid(cols = vars(code), rows = vars(Surface)) +
    scale_fill_gradient(low = "green", high = "red", name = "Value") +
    theme_minimal() +
    ylab("Tile") +
    xlab("Lane") +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    ggtitle("Cluster Density")
  
  
  ggplot(subset(tiles_df_long, code %in% c("Clusters", "Clusters_PF")), aes(x = as.factor(lane), y = as.factor(tile), fill = value)) +
    geom_tile() +
    facet_grid(cols = vars(code), rows = vars(Surface)) +
    scale_fill_gradient(low = "green", high = "red", name = "Value") +
    theme_minimal() +
    ylab("Tile") +
    xlab("Lane") +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    ggtitle("Clusters")
  
  
  ggplot(subset(tiles_df_long, code %in% c("PercAligned_R1", "PercAligned_R4")), aes(x = as.factor(lane), y = as.factor(tile), fill = value)) +
    geom_tile() +
    facet_grid(cols = vars(code), rows = vars(Surface)) +
    scale_fill_gradient(low = "green", high = "red", name = "Value") +
    theme_minimal() +
    ylab("Tile") +
    xlab("Lane") +
    theme(strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    ggtitle("Percent Aligned")
}




# Cluster Plot NextSeq
if (run_info$Sequencer == "NextSeq") {
  p1 <- tiles_df_long %>%
    filter(code == "Density" & Surface == "Top") %>%
    filter(lane == 1 | lane == 2) %>%
    mutate(tile = substring(tile, first = 4)) %>%
    ggplot(aes(x = Swath, y = as.factor(tile), fill = value / 10^3)) +
    geom_tile() +
    facet_grid(rows = vars(Camera), cols = vars(lane)) +
    theme_minimal() +
    ggtitle("Top Surface") +
    ylab("Tile") +
    labs(fill = "Density") +
    theme(legend.title = element_text(size = 12, face = "bold"),
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  p2 <- tiles_df_long %>%
    filter(code == "Density" & Surface == "Top") %>%
    filter(lane == 3 | lane == 4) %>%
    mutate(tile = substring(tile, first = 4)) %>%
    ggplot(aes(x = Swath, y = as.factor(tile), fill = value / 10^3)) +
    geom_tile() +
    facet_grid(rows = vars(Camera), cols = vars(lane)) +
    theme_minimal() +
    ggtitle("Top Surface") +
    ylab("Tile") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  p3 <- tiles_df_long %>%
    filter(code == "Density" & Surface == "Bottom") %>%
    filter(lane == 1 | lane == 2) %>%
    mutate(tile = substring(tile, first = 4)) %>%
    ggplot(aes(x = Swath, y = as.factor(tile), fill = value / 10^3)) +
    geom_tile() +
    facet_grid(rows = vars(Camera), cols = vars(lane)) +
    theme_minimal() +
    ggtitle("Bottom Surface") +
    ylab("Tile") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  p4 <- tiles_df_long %>%
    filter(code == "Density" & Surface == "Bottom") %>%
    filter(lane == 3 | lane == 4) %>%
    mutate(tile = substring(tile, first = 4)) %>%
    ggplot(aes(x = Swath, y = as.factor(tile), fill = value / 10^3)) +
    geom_tile() +
    facet_grid(rows = vars(Camera), cols = vars(lane)) +
    theme_minimal() +
    ggtitle("Bottom Surface") +
    ylab("Tile") +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "right")
}


# Clusters Boxplots
tiles_df_long %>%
  filter(code == "Clusters" | code == "Clusters_PF") %>%
  plot_ly(x = ~lane, y = ~value, color = ~code, type = "box") %>%
  layout(boxmode = "group",
         title = "Cluster",
         xaxis = list(title = "Lane", zeroline = FALSE, showgrid = TRUE), 
         yaxis = list(title = "Cluster Count", zeroline = FALSE, showgrid = TRUE))


# Density Boxplots
tiles_df_long %>%
  filter(code == "Density" | code == "Density_PF") %>%
  plot_ly(x = ~lane, y = ~value, color = ~code, type = "box") %>%
  layout(boxmode = "group",
         title = "Cluster Density",
         xaxis = list(title = "Lane", zeroline = FALSE, showgrid = TRUE), 
         yaxis = list(title = "Cluster Density", range = c(0, 250*10^3), zeroline = FALSE, showgrid = TRUE))
  

# Phasing Boxplots
tiles_df_long %>%
  filter(code %in% c("Phasing_R1", "Prephasing_R1", "Phasing_R4", "Prephasing_R4")) %>%
  plot_ly(x = ~lane, y = ~value*100, color = ~code, type = "box") %>%
  layout(boxmode = "group",
         title = "Phasing",
         xaxis = list(title = "Lane", zeroline = FALSE, showgrid = TRUE), 
         yaxis = list(title = "Phasing Rate %", zeroline = FALSE, showgrid = TRUE))


# Cleanup
rm(tiles_df, 
   tiles_df_long,
   p1, p2, p3, p4)





####### Old Code

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
