# Load necessary packages
library(DESeq2)
library(ggplot2)  # For plotting, if needed

# Define paths to files
original_metadata_path <- "/home/juliaa/genomanalys/Data/Metadata/Metadata.csv"
count_files_path <- "/home/juliaa/genomanalys/Data/Differential_expression/htseq/musang_king"

# Load metadata from the CSV file
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)
output_dir <- "/home/juliaa/genomanalys/Data/Differential_expression/output_musang_king"


# Load count files
count_files <- list.files(
    path = count_files_path,
    pattern = "*.txt",
    full.names = TRUE
)

# Check count files and their names
if (length(count_files) == 0) {
    stop("No count files found in the specified directory.")
}
print("Loaded count files:")
print(count_files)

# Read count data
count_data_list <- lapply(count_files, function(file) {
    # Read counts data from each file
    counts <- read.table(file, header = FALSE, col.names = c("gene", "counts"))
    
    # Filter out rows that do not represent typical gene identifiers
    valid_genes <- !grepl("^__", counts$gene)
    counts <- counts[valid_genes, ]
    
    # Convert counts to a named vector (gene counts)
    setNames(counts$counts, counts$gene)
})

# Extract base file names (sample identifiers) and assign to count data list
count_file_names <- gsub(".txt$", "", basename(count_files))
names(count_data_list) <- count_file_names

# Check the names of count files and their correspondence to metadata samples
print("Count file names (sample identifiers):")
print(count_file_names)
print("Metadata sample names:")
print(row.names(metadata))

# Check if there is any mismatch between count file names and metadata sample names
mismatch <- setdiff(row.names(metadata), count_file_names)
if (length(mismatch) > 0) {
    print("Mismatch between metadata and count files:")
    print(mismatch)
    # Optional: Handle mismatched samples in metadata, or count data as per your needs
    metadata <- metadata[row.names(metadata) %in% count_file_names,]
}

# Create count matrix
count_matrix <- do.call(cbind, count_data_list)

# Filter count matrix to include only samples in the filtered metadata
count_matrix <- count_matrix[, row.names(metadata), drop = FALSE]

# Convert tissue column to a factor (if necessary)
# Convert tissue column to a factor and set the reference level to "aril"
metadata$tissue <- factor(metadata$tissue, levels = c("Aril", "Root", "Leaf", "Stem"))

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ tissue)

# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
#res <- results(dds)

# Save the results to a CSV file
#write.csv(as.data.frame(res), file.path(output_dir, "01_DESeq2_results_mk.csv"))

# Optional: Plotting or additional analyses can be added here
# For example, creating an MA plot
#plotMA(dds, main = "MA-plot")

# Optional: Save the filtered metadata if needed
# write.csv(metadata, "/path/to/filtered_metadata.csv")

# Creating plots, such as MA plot
#plotMA(dds, main="MA-plot, treatment vs control")  # 'plotMA' function should take 'dds', not 'res'

# Specify the tissues to compare with aril
other_tissues <- c("Root", "Leaf", "Stem")

for (tissue in other_tissues) {
    # Specify the contrast: other tissue vs aril
    contrast <- c("tissue", tissue, "Aril")

    # Extract results for the specified contrast
    res <- results(dds, contrast = contrast)

    # Save the results to a CSV file
    write.csv(as.data.frame(res), file.path(output_dir, "02_DESeq2_results_mk.csv"))

    # Print some summary statistics
    print(paste("Results for comparison:", tissue, "vs aril"))
    print(summary(res))

    # Optional: Create an MA plot for the comparison
    plotMA(res, main = paste("MA-plot:", tissue, "vs aril"))
}


# Preparing data for the volcano plot
res_df <- as.data.frame(res)  # Convert 'res' to a data frame for plotting

# Generate the volcano plot using ggplot2
library(ggplot2)
# Dynamically set the upper limit based on the maximum value found, add a small buffer
max_padj <- max(-log10(res_df$padj), na.rm = TRUE)
min_padj <- min(-log10(res_df$padj), na.rm = TRUE)
upper_limit <- ceiling(max_padj + 1)

# Filter significant genes based on the given thresholds
significant_genes <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# Calculate the upper limit for the y-axis based on the maximum -log10(padj) value
max_padj <- max(-log10(res_df$padj), na.rm = TRUE)
upper_limit <- ceiling(max_padj + 1)

# Create a new column in res_df to categorize genes based on significance and log2 fold change direction
res_df$significance <- ifelse(
    res_df$padj < 0.05 & res_df$log2FoldChange > 1,
    "Upregulated",
    ifelse(
        res_df$padj < 0.05 & res_df$log2FoldChange < -1,
        "Downregulated",
        "Non-significant"
    )
)

# Define colors for each significance category
significance_colors <- c("Upregulated" = "red", "Downregulated" = "blue", "Non-significant" = "grey")

# Generate the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    # Plot the points using colors based on the significance column
    geom_point(alpha = 0.4) +
    # Apply the custom color palette
    scale_color_manual(values = significance_colors) +
    # Set x and y axis limits
    xlim(c(-8, 8)) +
    ylim(0, upper_limit) +
    # Add labels and title
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot comparing cultivars") +
    # Customize the plot appearance
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")

# Save the plot as a PNG file
output_path <- file.path(output_dir, "03_volcano_plot_mk.png")
ggsave(output_path, plot = volcano_plot, width = 10, height = 8)



library(pheatmap)

# Adjusted thresholds
new_padj_threshold <- 0.0000001  # Higher p-value threshold
new_log2FoldChange_threshold <- 6  # Higher log2FoldChange threshold

# Filter significant genes based on the new thresholds
significant_genes <- res_df[res_df$padj < new_padj_threshold & abs(res_df$log2FoldChange) > new_log2FoldChange_threshold, ]


# Filter out rows with NA values from significant_genes
significant_genes_filtered <- significant_genes[complete.cases(significant_genes), , drop = FALSE]

# Get row names of filtered significant genes
significant_genes_list <- rownames(significant_genes_filtered)


# Extract normalized counts from the DESeqDataSet
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(normalized_counts, file.path(output_dir, "01norm_counts_aril.csv"))


# Subset the normalized counts matrix to include only the filtered significant genes
filtered_normalized_counts <- normalized_counts[significant_genes_list, , drop = FALSE]

write.csv(filtered_normalized_counts, file.path(output_dir, "01norm_counts_filtered_aril.csv"))

sample_names <- colnames(filtered_normalized_counts)
tissue_info <- metadata$tissue[match(sample_names, rownames(metadata))]
colnames(filtered_normalized_counts) <- paste(sample_names, tissue_info, sep = "_")


# Generate the heatmap
heatmap_plot <- pheatmap(filtered_normalized_counts,
                         scale = "row",  # Scale rows (genes) for better visualization
                         clustering_method = "average",  # Clustering method (e.g., "average", "ward.D2")
                         cluster_cols = TRUE,  # Cluster samples (columns)
                         cluster_rows = TRUE,  # Cluster genes (rows)
                         show_rownames = TRUE,  # Show gene names (row names)
                         show_colnames = TRUE,  # Show sample names (column names)
                         fontsize_row = 5,  # Font size for gene names
                         fontsize_col = 8,  # Font size for sample names
                         main = "Heatmap of Normalized Counts")  # Title of the heatmap

# Save the heatmap as a PNG file
output_heatmap_path <- file.path(output_dir, "06_heatmap_mk.png")
png(output_heatmap_path, width = 10, height = 8, units = "in", res = 300)
print(heatmap_plot)
dev.off()

