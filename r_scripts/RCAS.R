if (!requireNamespace("BiocManager", quietly=TRUE)) 
  install.packages("BiocManager")
BiocManager::install('RCAS')

install.packages("ggplot2")

library(RCAS)
library(ggplot2)

queryRegions <- importBed(filePath = "initial_peaks.bed", sampleN = 10000)

# Annotationfile: hg38_UCSC.gtf (downloaded from UCSC Genome Browser)
gff <- importGtf(filePath = "hg38_UCSC.gtf")

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
ggsave(fig,filename = "RCAS_result_UCSC.png",width = 10,height = 8)

# Annotation file: hg38_ENSEMBL.gtf (downloaded from ENSEMBLE database)
gff <- importGtf(filePath = "hg38_ENSEMBL.gtf")

overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))

biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
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
ggsave(fig,filename = "RCAS_result_ENSEMBL.png",width = 10,height = 8)
