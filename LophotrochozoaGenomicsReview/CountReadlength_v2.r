library(dplyr)
library(ggplot2)
library(introdataviz)

Rawreads <- read.table("~/000_changed_documents_20250823/Analyses/InvertOmics/Gnathostomulida/Gnathostomula_New/Gnathi_2runs.fastq_CountReadlength.txt", quote="\"", comment.char="")

Rawreads$X1 <- 'all_reads'

# N50 calculation function
calculate_N50 <- function(lengths) {
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total <- sum(sorted_lengths)
  cum_sum <- cumsum(sorted_lengths)
  n50 <- sorted_lengths[which(cum_sum >= total * 0.5)[1]]
  return(n50)
}

read_stats <- function(y, upper_limit = max(Rawreads$V1) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "Count =", length(y), "\n",
      "Mean =", round(mean(y), 0), "\n",
      "Median =", round(median(y), 0), "\n",
      "N50 =", calculate_N50(y), "\n",
      "Amount (in Gb) =", round(sum(y)/1000000000, 2), "\n"
    )
  ))
}

ggplot(Rawreads, aes(x=X1, y=V1, fill = X1)) + 
  introdataviz::geom_split_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  stat_summary(fun.data = read_stats, geom = "text", hjust = 0, vjust = 0.9, position = position_dodge(0))

