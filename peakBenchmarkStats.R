library(DescTools)
library(dplyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(gridExtra)

measures <- read.csv("measures.txt", sep = "\t")
names(measures)[1] <- "Peak.Caller"

clean_df <- function(df, name){
  dd <- df[grepl(paste0(name, "[.]"), df[[1]]),]
  dd[[1]] <- gsub("[-].*$", "", dd[[1]])
  dd[[1]] <- gsub(paste0("^", name, "[.]"), "", dd[[1]])
  return(dd)
}

find_cutoff <- function(df){
  df$cutoff <- gsub("^.*[.].*[-]", "", df[[1]])
  df$cutoff <- gsub("_peaks", "", df$cutoff)
  return(df)
}

get_AUC <- function(z){
  df <- data.frame(Doubles=double())
  names(df)[1] <- as.character(substitute(z))
  for(val in unique(z[[1]])){
    AUC <- AUC(x=z[ which(z$Peak.Caller == val), ][['Recall']],
               y=z[ which(z$Peak.Caller == val), ][['Precision']],
               method = 'trapezoid')
    df[paste0(val),] <- AUC
  }
  return(df)
}

Omni <- clean_df(measures, "Omni")
Original <- clean_df(measures, "Original")
Parker <- clean_df(measures, "Parker")

plotPR_Curve <- function(name){
  ggplot(name, aes(x=Recall, y=Precision, group=Peak.Caller)) +
    geom_line(aes(color=Peak.Caller))+
    geom_point(aes(color=Peak.Caller))+
    scale_color_discrete(name = "Peak Calling Method",
                         labels = c("BAM", "BAMPE", "BED", "Genrich",
                                    "Genrich (Replicate Mode)",
                                    "HMMRATAC (summits += 50bp)",
                                    "HMMRATAC (gappedPeak"))
}

prow <- plot_grid(
  plotPR_Curve(Omni) + theme(legend.position="none"),
  plotPR_Curve(Original) + theme(legend.position="none"),
  plotPR_Curve(Parker) + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  plotPR_Curve(Omni) + theme(legend.box.margin = margin(0, 0, 0, 0))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot::plot_grid(prow, legend, rel_widths = c(2, .4))

AUC <- cbind(get_AUC(Omni), get_AUC(Original), get_AUC(Parker))
table <- grid.table(AUC, rows=c( "BAM", "BAMPE", "BED", "Genrich",
                               "Genrich.Replicates", "HMMRATAC.gapped_peaks",
                               "HMMRATAC.summits"))


#====Coverage Statistics====
coverageStatistics <- read.csv("coverageStatistics_sorted.txt", sep = "\t")
names(coverageStatistics)[1] <- "Peak.Caller"

tidyCoverage <- function(df, a, b){
  df <- find_cutoff(df)
  df$Tool <- Parker$Peak.Caller
  df <- reshape::melt(df, id=c("Peak.Caller", "Tool", "cutoff"))%>%
    dplyr::filter(cutoff %in% c(a, b))# %>%
    #dplyr::filter(!((Tool %in% c("HMMRATAC", "HMMRATAC_peaks")) &
    #                  (cutoff == "p05"))) %>%
   # dplyr::filter(!((Tool %in% c("Genrich", "Genrich.Replicates")) &
    #                  (cutoff == "p05")))
  return(df)
}

plotCoverage <- function(df){
  labels = c("BAM", "BAMPE", "BED", "Genrich",
             "Genrich (Replicate Mode)",
             "HMMRATAC (summits += 50bp)",
             "HMMRATAC (gappedPeak")
  ggplot(data = df, aes(x=Tool, y=value, fill=variable)) +
    geom_bar(stat="identity") +
    xlab("Peak Calling Method") +
    ylab("Number of Base Pairs") +
    labs(fill = "Coverage") +
    scale_fill_hue(labels = c("0-2 bp", "3-6 bp", ">6 bp")) +
    scale_x_discrete(labels = labels)
}

Omni.cov <- tidyCoverage(coverageStatistics[1:55,], "p01", "p30")
Original.cov <- tidyCoverage(coverageStatistics[152:206,], "p01", "p30")
Parker.cov <- tidyCoverage(coverageStatistics[351:405,], "p01", "p30")

plotCoverage(Omni.cov)
plotCoverage(Original.cov)
plotCoverage(Parker.cov)
plotCoverage(tidyCoverage(coverageStatistics[1:55,], "p1", "p20"))
plotCoverage(tidyCoverage(coverageStatistics[152:206,], "p1", "p20"))
plotCoverage(tidyCoverage(coverageStatistics[351:405,], "p1", "p20"))

