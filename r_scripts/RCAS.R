if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install('RCAS')
install.packages("ggplot2")
install.packages("getopt")

library(RCAS)
library(ggplot2)
library(getopt)

spec <- matrix(
  c("input_bed",  "b", 2, "character",
    "input_gtf", "g", 2, "character",  
    "output_plot", "o", 2, "character", 
    "help",   "h", 2, "logical"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec=spec)

input_bed = opt$input_bed
input_gtf = opt$input_gtf
output_plot = opt$output_plot

# get input .bed file
queryRegions <- importBed(filePath = input_bed, sampleN = 10000)
# get input .gtf file
gff <- importGtf(filePath = input_gtf)

overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))

biotype_col <- grep('type', colnames(overlaps), value = T)
df <- overlaps[,length(unique(queryIndex)), by = biotype_col]
colnames(df) <- c("feature", "count")
df$percent <- round(df$count / length(queryRegions) * 100, 1)
df <- df[order(count, decreasing = TRUE)]
fig = ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 0.5), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))

# output the plot
ggsave(fig,filename = output_plot,width = 10,height = 8)

