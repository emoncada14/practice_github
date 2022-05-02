install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")

#Load the packages

library(qiime2R)

library(tidyverse)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/core-metrics-results/* core-metrics-results/.
##############################################


setwd("~/Desktop/Molecular Microbiome/Final Project/New_final_project/")
#setwd("/Users/john2185/Documents/")
list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("finalmetadataEM2.tsv")
#metadata2 <- read.delim("sample-metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#metadata2 <- metadata2[-1,]
str(metadata)
levels(metadata$`treatment`)


row.names(metadata) <- metadata[ ,1]
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")

#gut.region_colors <- c("Black", "Blue")
treatment_colors <- c("Black", "Blue", "Green", "Gray")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

my_column <- "gut.region"
my_column <- "treatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=gut.region_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a tiff 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~treatment)
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=gut.region_colors, name = my_column)
ggsave(paste0("output/BC-basic_gutregionbygroups", my_column,".pdf"), height=3, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a tiff 3 inches by 4 inches


centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "gut.region"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=treatment_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_treatment", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,"-gut.regionPP.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

##SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

#gut.region_colors <- c("Black", "Blue")
treatment_colors <- c("Black", "Blue", "Green", "Gray")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

#my_column <- "gut.region"
my_column <- "treatment"

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=gut.region_colors, name = "Gut Region")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=gut.region_colors, name = "Gut Region")
ggsave(paste0("output/Wuni-ellipse_", my_column,"trt.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Jaccard

Jaccard_PCoA<-read_qza("core-metrics-results/jaccard_pcoa_results.qza")

#gut.region_colors <- c("Black", "Blue")
treatment_colors <- c("Black", "Blue", "Green", "Gray")

jaccard_meta <- Jaccard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

#my_column <- "gut.region"
my_column <- "treatment"

ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=treatment_colors, name = my_column)
ggsave(paste0("output/Jaccard-", my_column,".pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jaccard_meta,mean)
colnames(centroids)[1] <- "treatment"

ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=treatment_colors, name = my_column)
ggsave(paste0("output/Jaccard-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= treatment)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("output/Jaccard-ellipse_", my_column,"-gut.regionPP.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

#Unweighted Unifrac
UWuni_PCoA<-read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

#gut.region_colors <- c("Black", "Blue")
treatment_colors <- c("Black", "Blue", "Green", "Gray")

UWuni_meta <- UWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

#my_column <- "gut.region"
my_column <- "treatment"

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UWuni_meta,mean)

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=treatment_colors, name = "Treatment")
ggsave(paste0("output/UWUni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= gut.region), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=treatment_colors, name = "Treatment")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


© 2022 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
