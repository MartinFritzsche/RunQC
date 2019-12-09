####### Quality Metric ######

# Master: RunQC_script.R

# Parse quality metrics (depends on sequencer)
if (run_info$Sequencer == "MiSeq") {
  quality_df <- qualityMetrics(fc) # This does not work for NextSeq!
} else {
  quality_df <- fc@parsedData$savQualityFormatV5@data # This does!
}


###### QScore Distribution Overview ######
quality_overview <- quality_df %>%
  select(Q1:Q47) %>%
  colSums() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Quality") %>%
  mutate(Quality = as.integer(substring(Quality, 2))) %>%
  mutate(Threshold = ifelse(Quality > 29, "Pass", "Fail")) %>%
  rename(Total = '.') %>%
  mutate(Perc = Total / sum(Total))


# Q Score Distribution
ggplot(quality_overview, aes(x = as.factor(Quality), y = Total / 10^9, fill = Threshold)) +
  geom_bar(stat = "identity") +
  xlab("Q Score") +
  ylab("Total Gb") +
  ggtitle("Q Score Distribution") +
  theme_minimal() +
  theme(legend.position = "none")


# Same with plotly
quality_overview %>%
  plot_ly(x = ~as.factor(Quality), y = ~Total, color = ~Threshold, 
          hoverinfo = "text", text = ~paste0(round(Perc * 100, 2), "%"),
          colors = c("firebrick2", "darkolivegreen2")) %>%
  add_bars() %>%
  layout(title = "Q Score Distribution",
         showlegend = FALSE, 
         xaxis = list(title = "Q Score", dtick = 1, tickangle = 0), 
         yaxis = list(title = "Total Bases"))

rm(quality_overview)

###### Quality overview by cycle ######

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


# Same with plotly
quality_overview_cycle %>%
  plot_ly(x = ~cycle, y = ~Perc_Pass,
          hoverinfo = "text", 
          text = ~paste("Tile: ", tile, "<br>",
                        "Cycle: ", cycle)) %>%
  add_boxplot(boxpoints = "outliers", 
              marker = list(color = "red")) %>%
  layout(title = "Quality Overview by Cycle",
         xaxis = list(title = "Cycle", dtick = 10, tickangle = 90), 
         yaxis = list(title = "% Clusters >= Q30"))


###### Quality overview by tile and cycle ######

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


