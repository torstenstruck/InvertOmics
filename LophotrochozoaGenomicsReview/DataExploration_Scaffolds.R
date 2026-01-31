#Analysis of genome data of Lophotrochozoa

library(ggplot2) #generating violin plots
library(dplyr) #calculating statistics
library(tidyr) #generating pivot tables
library(writexl) #write excel file
library(readxl) #read excel file
library(tidyverse)

# Set working environment
setwd("~/000_changed_documents_20250823/Analyses/InvertOmics/01_NCBI_GoaT_data_20251203")
dir.create("ResultsScaffolds")

#import genome data from GoaT 
Species_Data_Input <- read_excel("./LophotrochozoaGenomesScaffoldsGoaT.xlsx")
Species_Data_Input <- as.data.frame(Species_Data_Input)

# Count occurrences of phylum and conditional counts for annotation file available as well
phylum_counts <- Species_Data_Input %>%
  group_by(phylum) %>%
  summarise(
    all = n(),  # Total count of each phylum
    annotated = sum(!is.na(gene_count) & gene_count != 0)  # Count only if count of genes has valid values
  ) %>%
  arrange(desc(all))

# Reshape the data for plotting
phylum_counts_long <- phylum_counts %>%
  pivot_longer(cols = c(all, annotated), names_to = "Type", values_to = "Value")

# Order the Phylum factor by the 'all' counts after reshaping
phylum_counts_long$phylum <- factor(phylum_counts_long$phylum, levels = phylum_counts$phylum[order(-phylum_counts$all)])

# Create a dodged bar plot with both counts displayed next to each other
phylum_plot <- ggplot(phylum_counts_long, aes(x = phylum, y = Value, fill = Type)) +  
  geom_bar(stat = "identity", position = position_dodge()) +  
  geom_text(aes(label = Value), position = position_dodge(0.9), vjust = -0.5, color = "black") +  
  labs(title = "Count of Genomes per Phylum", x = "Phylum", y = "Number of Genomes") +
  scale_fill_manual(values = c("skyblue", "orange"), labels = c("Total Count", "Count Annotated")) +  
  theme_minimal()

# Print phylum-level plot
print(phylum_plot)

# Save the phylum-level plot
ggsave(filename = paste0("ResultsScaffolds/Counts_by_Phylum.pdf"), 
       plot = phylum_plot, 
       device = "pdf", 
       width = 29.7,  # Width for A4 in cm
       height = 21,   # Height for A4 in cm
       units = "cm")  # Set units to centimeters

#### Analyse the numbers within the phyla
# Initialize a list to store results and a data frame for final results
results_list <- list()
final_results <- data.frame()

# Loop through each unique Phylum and create separate results
unique_phyla <- unique(Species_Data_Input$phylum)

for (Phylum in unique_phyla) {
  # Filter the data for current phylum
  df <- Species_Data_Input %>% filter(phylum == Phylum)
  
  # Class-level analysis: count total and annotated occurrences
  class_counts <- df %>%
    group_by(class) %>%
    summarise(
      all = n(),  # Total occurrences of each Class
      annotated = sum(!is.na(gene_count) & gene_count != 0)  # Count only if count of genes has valid values
    )
  
  # Add phylum name for clarity
  class_counts <- class_counts %>%
    mutate(phylum = Phylum)
  
  # Store the results in the results_list
  results_list[[gsub(" ", "_", Phylum)]] <- class_counts
  
  # Reshape for plotting
  class_counts_long <- class_counts %>%
    pivot_longer(cols = c(all, annotated), names_to = "Type", values_to = "Count")
  
  # Order the Class factor by the 'all' counts after reshaping
  class_counts_long$class <- factor(class_counts_long$class, levels = class_counts$class[order(-class_counts$all)])
  
  # Create a dodge bar plot for Class counts for the current phylum
  p <- ggplot(class_counts_long, aes(x = class, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = Count), position = position_dodge(0.9), vjust = -0.5) +
    labs(title = paste("Counts by Class for ", Phylum), 
         x = "Class", y = "Count") +
    theme_minimal() +
    scale_fill_manual(values = c("skyblue", "orange"), labels = c("Total Count", "Count Annotated"))
  
  # Print the plot to ensure it is displayed
  print(p) 
  
  # Save each plot as a PDF in A4 landscape format
  ggsave(filename = paste0("ResultsScaffolds/Counts_by_Class_", gsub(" ", "_", Phylum), ".pdf"), 
         plot = p, 
         device = "pdf", 
         width = 29.7,  # Width for A4 in cm
         height = 21,   # Height for A4 in cm
         units = "cm")  # Set units to centimeters
  
  # Combine results into the final_results data frame
  final_results <- bind_rows(final_results, class_counts)
}

# Adjust the order of the columns for better presentation
final_results <- final_results[, c("phylum", "class", "all", "annotated")]

# Write the final results to an Excel file
write_xlsx(final_results, path = "ResultsScaffolds/Class_Counts_By_Phylum.xlsx")

# Print final results
print(final_results)

### Analyze the distribution of genomes across genera and species
# Function to analyze a specific column
analyze_distribution <- function(data, column_name) {
  # Count occurrences of each unique value
  occurrence_counts <- data %>%
    group_by(!!sym(column_name)) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count))
  
  # Include Phylum information for the first occurrences
  phylum_info <- data %>%
    group_by(!!sym(column_name)) %>%
    summarise(Phylum = first(phylum)) %>%
    ungroup()
  
  # Merge phylum info with occurrence_counts
  occurrence_counts <- left_join(occurrence_counts, phylum_info, by = column_name)
  
  # Rearrange columns into different order
  occurrence_counts <- occurrence_counts %>%
    select(Phylum, !!sym(column_name), Count)
  print(occurrence_counts)
  
  # Export to Excel
  write_xlsx(occurrence_counts, path = paste0("ResultsScaffolds/", column_name, "_occurrences.xlsx"))
  
  # Count how many genera/species have the same occurrence count
  occurrence_distribution <- occurrence_counts %>%
    group_by(Count) %>%
    summarise(Occurrence_Frequency = n()) %>%
    arrange(desc(Count))
  
  # Generate the bar plot for occurrence distribution
  plot <- ggplot(occurrence_distribution, aes(x = Count, y = Occurrence_Frequency)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Occurrence_Frequency), vjust = -0.5) +  # Add labels above the bars
    labs(title = paste("Distribution of Occurrences of", column_name),
         x = "Occurrences",
         y = paste0("Number of ", column_name)) +
    theme_minimal()
  
  # Print plot
  print(plot)
  
  # Save plot
  ggsave(filename = paste0("ResultsScaffolds/", column_name, "_occurrences_plot.pdf"),
         plot = plot, 
         device = "pdf", 
         width = 29.7,  # Width for A4 in cm
         height = 21,   # Height for A4 in cm
         units = "cm")  # Set units to centimeters
}

# Analyze Genus
analyze_distribution(Species_Data_Input, "genus")
# Analyze Species
analyze_distribution(Species_Data_Input, "species")

### Analyse all parameters for each assembly
# Loop through each relevant column
for (i in 8:29) {
  # Access the column using its index
  column <- colnames(Species_Data_Input)[i]
  
  # Calculate max and min for log scale decision
  column_min <- min(Species_Data_Input[[column]], na.rm = TRUE)
  column_max <- max(Species_Data_Input[[column]], na.rm = TRUE)
  log_scale <- (log10(column_max) - log10(column_min)) > 1  # Checks for more than two order of magnitude
  
  ### Plot for All Data ###
  # Create ggplot for overall data
  p_all_data <- ggplot(Species_Data_Input, aes_string(x = "1", y = column)) + # "1" to treat all data as one group
    geom_boxplot(alpha = 0.5) +  # Boxplot for all data
    geom_violin(alpha = 0.3) +  # Violin plot for all data
    labs(title = paste("Box and Violin Plot for All Data -", column),
         x = "All",
         y = column) +
    theme_minimal()
  
  # If log scale is applicable
  if (log_scale) {
    p_all_data <- p_all_data + scale_y_log10()  # Set y-axis to logarithmic scale
  }
  
  # Calculate summary statistics for overall data
  summary_stats_all <- data.frame(
    Statistic = c("Min", "Q1", "Median", "Q3", "Max"),
    Value = c(column_min, 
              quantile(Species_Data_Input[[column]], 0.25, na.rm = TRUE, type = 1), 
              median(Species_Data_Input[[column]], na.rm = TRUE), 
              quantile(Species_Data_Input[[column]], 0.75, na.rm = TRUE, type = 1), 
              column_max)
  )
  
  # Add summary stats as text annotations
  p_all_data <- p_all_data +
    geom_text(data = summary_stats_all, aes(x = 1, y = Value, label = Value), hjust = -0.5, vjust = 0)
  
  # Save the overall data plot
  ggsave(filename = paste0("ResultsScaffolds/", column, "_All_Data_Plot.pdf"), plot = p_all_data, 
         device = "pdf", 
         width = 10,  # Width for A4/3 in cm
         height = 21,   # Height for A4 in cm
         units = "cm")  # Set units to centimeters
  
  ### Calculate and Save Summary Statistics for Grouped Data by Phylum ###
  summary_stats_grouped <- Species_Data_Input %>%
    group_by(phylum) %>%
    summarise(
      Count = sum(!is.na(get(column))),  # Count of non-NA entries in the specified column
      Max = max(get(column), na.rm = TRUE),
      Q1 = quantile(get(column), 0.25, na.rm = TRUE, type = 1),
      Median = median(get(column), na.rm = TRUE),
      Q3 = quantile(get(column), 0.75, na.rm = TRUE, type = 1),
      Min = min(get(column), na.rm = TRUE)
    )
  
  # Save the summary statistics to a separate Excel file for each column
  write_xlsx(summary_stats_grouped, path = paste0("ResultsScaffolds/", column, "_Summary_Statistics_By_Phylum.xlsx"))
  
  # Filter to keep only phyla with at least 5 entries
  plot_data <- summary_stats_grouped %>% filter(Count >= 5)
  
  ### Plot for Grouped by Phylum ###
  # Check if there are any remaining entries for plotting
  if (nrow(plot_data) > 0) {
    # Create ggplot for grouped by Phylum
    p_grouped_phylum <- ggplot(Species_Data_Input[Species_Data_Input$phylum %in% plot_data$phylum, ], 
                               aes_string(x = "phylum", y = column)) +
      geom_boxplot(aes(fill = phylum), alpha = 0.5) +  # Boxplot grouped by Phylum
      geom_violin(aes(fill = phylum), alpha = 0.3) +  # Violin plot grouped by Phylum
      labs(title = paste("Box and Violin Plot Grouped by Phylum -", column),
           x = "Phylum",
           y = column) +
      theme_minimal()
    
    # If log scale is applicable
    if (log_scale) {
      p_grouped_phylum <- p_grouped_phylum + scale_y_log10()  # Set y-axis to logarithmic scale
    }
    
    # Add summary stats as text annotations for each Phylum
    p_grouped_phylum <- p_grouped_phylum +
      geom_text(data = plot_data, aes(x = phylum, y = Median, label = round(Median, 2)), hjust = -0.3, vjust = 0) +
      geom_text(data = plot_data, aes(x = phylum, y = Q1, label = round(Q1, 2)), hjust = -0.3, vjust = 0) +
      geom_text(data = plot_data, aes(x = phylum, y = Q3, label = round(Q3, 2)), hjust = -0.3, vjust = 0) +
      geom_text(data = plot_data, aes(x = phylum, y = Max, label = round(Max, 2)), hjust = -0.3, vjust = 0, color = "blue") +
      geom_text(data = plot_data, aes(x = phylum, y = Min, label = round(Min, 2)), hjust = -0.3, vjust = 0, color = "blue") 
    
    # Save the grouped by Phylum plot
    ggsave(filename = paste0("ResultsScaffolds/", column, "_Grouped_by_Phylum_Plot.pdf"), plot = p_grouped_phylum, 
           device = "pdf", 
           width = 29.7,  # Width for A4 in cm
           height = 21,   # Height for A4 in cm
           units = "cm")  # Set units to centimeters
    
    # Show plots
    print(p_grouped_phylum)
  }
  
  # Show plots
  print(p_all_data)
}

#Analyse different parameters against each other and color by phylum
analyze_linear_regression <- function(raw_data, color, column_x, column_y) {
  data <- cbind(raw_data[[color]],raw_data[[column_x]],raw_data[[column_y]])
  data <- na.omit(data)
  data <- as.data.frame(data)
  data$V2 <- as.numeric(data$V2)
  data$V3 <- as.numeric(data$V3)

  #Fit a linear model
  model_est_ass <- lm(data$V3 ~ data$V2)

  #Extract coefficients and R^2 value
  coef_model_est_ass <- coef(model_est_ass)
  r_squared_model_est_ass <- summary(model_est_ass)$r.squared

  # Checks for more than three orders of magnitude in both dimensions
  columnx_min <- min(data$V2, na.rm = TRUE)
  columnx_max <- max(data$V2, na.rm = TRUE)
  log_scalex <- (log10(columnx_max) - log10(columnx_min)) > 3  
  columny_min <- min(data$V3, na.rm = TRUE)
  columny_max <- max(data$V3, na.rm = TRUE)
  log_scaley <- (log10(columny_max) - log10(columny_min)) > 3  
  
  #Plot the data with trendline and confidence interval
  linreg_p <- ggplot(data, aes(x = V2, y = V3, colour = V1)) +
    geom_point() + # Scatter plot of the data
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +  # Trendline with CI
    geom_abline(col = "red", lty = 2) +  # Line where x = y
    labs(x = column_x, y = column_y) +
    annotate("text", x = 1, y = 10, label = paste("y = ", round(coef_model_est_ass[2], 2), "x + ", round(coef_model_est_ass[1], 2), "\nR^2 = ", round(r_squared_model_est_ass, 2)), 
         color = "black", size = 3, hjust = 0)  # Add formula and R2
  # If log scale is applicable
  if (log_scalex) {
    linreg_p <- linreg_p + scale_x_log10()  # Set x-axis to logarithmic scale
  }
  if (log_scaley) {
    linreg_p <- linreg_p + scale_y_log10()  # Set y-axis to logarithmic scale
  }
  
  #Distribution of the percentage deviation from the estimated genome size
  perdev <- ((data$V3 - data$V2)/data$V2*100)
  # Plot the distribution using a histogram
  hist_p <- hist(perdev, xlab=paste("Percentage deviation of ", column_y, " from ", column_x), col="skyblue", border="black", freq=FALSE, breaks=31)
  # Add a density curve
  hist_p <- lines(density(na.omit(perdev)), col="red", lwd=2)
  #Show plot
  print(hist_p)
  
  dev.copy(pdf, file = paste0("ResultsScaffolds/", column_y, "_vs_", column_x, "_Percentage_deviation.pdf"))
  dev.off()  

  #Show plot
  print(linreg_p)
  
  # Save the Percentage distribution plot
  ggsave(filename = paste0("ResultsScaffolds/", column_y, "_vs_", column_x, "_Linear_regression.pdf"), plot = linreg_p, 
         device = "pdf", 
         width = 29.7,  # Width for A4 in cm
         height = 21,   # Height for A4 in cm
         units = "cm")  # Set units to centimeters
  
}

# Analyze assembly_span vs. genome_size_estimate
analyze_linear_regression(Species_Data_Input, "phylum", "genome_size_estimate", "assembly_span")

# Analyze assembly_span vs. genome_size_estimate
analyze_linear_regression(Species_Data_Input, "phylum", "genome_size_direct", "assembly_span")
