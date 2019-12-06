input_demultiplex <- readLines(paste0(path, "/DemultiplexSummaryF1L1.txt"), n = 30)

demultiplex <- read.table(text = input_demultiplex, sep = "\t", skip = 1, header = TRUE)

demultiplex_tidy <- demultiplex %>%
  select(-c("X")) %>%
  rename(Tile = SampleName) %>%
  pivot_longer(cols = -c("Tile"), names_to = "Sample", values_to = "Percent")

demultiplex_tidy <- mutate(demultiplex_tidy, Tile = as.integer(substring(demultiplex_tidy$Tile, 5)))

ggplot(demultiplex_tidy, aes(x = Sample, y = Percent)) +
  geom_boxplot(colour = "blue", outlier.colour = "red") +
  coord_flip() +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = -0.25),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line()) +
  ylab("Percent Identified") +
  ggtitle("Sample Demultiplex")


ggplot(demultiplex_tidy, aes(x = as.factor(Tile), y = Sample, fill = Percent)) +
  geom_tile(colour = "gray33") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Tile") +
  ggtitle("Sample Demultiplex per Tile")





input_indexes <- readLines(paste0(path, "/DemultiplexSummaryF1L1.txt"))
input_indexes1 <- input_indexes[35:134]
input_indexes2 <- input_indexes[136:length(input_indexes)]

indexes1 <- read.table(text = input_indexes1, sep = "\t", header = FALSE)
indexes2 <- read.table(text = input_indexes2, sep = "\t", header = FALSE)

colnames(indexes1) <- c("Index", "RevComp_Index", "Count")
colnames(indexes2) <- c("Index", "RevComp_Index", "Count")

samples_sheet <- read.csv(paste0(path, "/SampleSheet.csv"), header = TRUE, skip = 17)

indexes1$InSampleSheet <- FALSE
indexes1$InSampleSheet[which(indexes1$Index %in% samples_sheet$index)] <- TRUE

indexes2$InSampleSheet <- FALSE
indexes2$InSampleSheet[which(indexes2$Index %in% samples_sheet$index2)] <- TRUE

ggplot(indexes1, aes(x = Index, y = Count, fill = InSampleSheet)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = Count), size = 2, hjust = -0.4) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = -0.08),
        axis.text.y = element_text(size = 10,  family = "Courier New"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top") +
  labs(fill = "Included in Sample Sheet?") +
  ggtitle("Top 100 Identified Index 1")

ggplot(indexes2, aes(x = Index, y = Count, fill = InSampleSheet)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(label = Count), size = 2, hjust = -0.4) +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = -0.08),
        axis.text.y = element_text(size = 10,  family = "Courier New"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top") +
  labs(fill = "Included in Sample Sheet?") +
  ggtitle("Top 100 Identified Index 2")


# which(indexes1$Index %in% samples_sheet$index)
# which(indexes1$RevComp_Index %in% samples_sheet$index)
# which(indexes1$Index %in% samples_sheet$index2)
# which(indexes1$RevComp_Index %in% samples_sheet$index2)
# 
# which(indexes2$Index %in% samples_sheet$index2)
# which(indexes2$RevComp_Index %in% samples_sheet$index2)
# which(indexes2$Index %in% samples_sheet$index1)
# which(indexes2$RevComp_Index %in% samples_sheet$index1)


## Test for Uniformity (goodness of fit) with Chi-Square

# test_uniformity_chi2 <- function(input, significance = 0.05) {
#   result = (chisq.test(input)$p.value > significance)
#   return(result)
# }
# 
# demultiplex_unitest <- demultiplex %>%
#   select(-c("SampleName", "X"))
# 
# demultiplex_unitest$Test <- rnorm(n = nrow(demultiplex_unitest), mean = 10, sd = 2)
# 
# demultiplex_unitest_result <- apply(demultiplex_unitest, 2, test_uniformity_chi2)
# 
# which(!demultiplex_unitest_result)


