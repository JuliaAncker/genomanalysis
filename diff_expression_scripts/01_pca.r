library(Deseq2)
library(ggplot2)  # For plotting, if needed

# Define paths to files
original_metadata_path <- "/home/juliaa/genomanalys/Data/Metadata/Metadata.csv"
count_files_path <- "/home/juliaa/genomanalys/Data/Differential_expression/htseq"

# Load metadata from the CSV file
metadata <- read.csv(original_metadata_path, header = TRUE, row.names = 1)
output_dir <- "/home/juliaa/genomanalys/Data/Differential_expression"


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
metadata$Cultivar <- as.factor(metadata$Cultivar)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ Cultivar)

# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

res_df <- as.data.frame(res)
colnames(res_df)[1] <- "gene"

# Save the results to a CSV file
write.csv(res_df, file.path(output_dir, "03_DESeq2_results_aril.csv"))



#PCA

# Normalize the counts data using variance stabilizing transformation (VST)
rlog_counts <- assay(rlog(dds, blind = FALSE))

# Run PCA using prcomp() on the transposed VST counts data
pca_result <- prcomp(t(rlog_counts))

# Extract the first two principal components
pc1 <- pca_result$x[, 1]
pc2 <- pca_result$x[, 2]

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pc1, PC2 = pc2, Sample = colnames(rlog_counts))

# Add metadata to the data frame for plotting
pca_df <- cbind(pca_df, metadata)

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Calculate the variance explained by PC1 and PC2
variance_PC1 <- round(variance_explained[1] * 100, 2)
variance_PC2 <- round(variance_explained[2] * 100, 2)


# Plot the PCA using ggplot2
# PCA Plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    # Add points colored based on the cultivar information in the metadata
    geom_point(aes(color = colData(dds)$tissue, shape = colData(dds)$Cultivar), size = 3) +
    # Add sample names as labels to the points using geom_text()
    geom_text(aes(label = colnames(counts(dds))), vjust = -0.5, hjust = 0.5, size = 3, color = "black") +
    # Add labels for the axes and title
    labs(title = "PCA Plot", x = paste("PC1 (", variance_PC1, "% variance)", sep = ""), 
         y = paste("PC2 (", variance_PC2, "% variance)", sep = ""),
         color = "Tissues",
         shape = "Cultivar") +
    expand_limits(x = c(min(pca_df$PC1) - 2, max(pca_df$PC1) + 2)) +
    # Use a minimalistic theme
    theme_bw()

output_path <- file.path(output_dir, "01_pca_plot.png")
ggsave(output_path, plot = pca_plot, width = 8, height = 6)
	
# Save and display
