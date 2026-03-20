library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(filename, delimiter) {
    df <- read.delim(filename, sep = delimiter, header = TRUE,
                     row.names = 1, check.names = FALSE)
    return(as.data.frame(df))
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
    # sdev holds the standard deviation of each PC;
    # squaring gives variance, dividing by the sum gives proportion
    var_explained <- (pca_results$sdev ^ 2) / sum(pca_results$sdev ^ 2)
    return(var_explained)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
    tibble(
        principal_components= colnames(pca_results$x),  
        variance_explained= pca_ve,
        cumulative= cumsum(pca_ve)
    )
}

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
    meta_df <- read.csv(metadata, row.names = 1)
 
    scores <- as.data.frame(pca_results$x)
    scores$SixSubTypesClassification <-
        meta_df[rownames(scores), "SixSubTypesClassification"]
 
    ggplot(scores, aes(x = PC1, y = PC2,
                       color = SixSubTypesClassification)) +
        geom_point(size = 2, alpha = 0.8) +
        labs(
            title = "PCA Biplot: PC1 vs PC2",
            x     = "PC1",
            y     = "PC2",
            color = "Six Subtype\nClassification"
        ) +
        theme_bw() +
        theme(legend.position = "right")
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
    sig <- diff_exp_tibble %>%
        filter(padj < fdr_threshold) %>%
        pull(probeid)           # first column holds the probe IDs
    return(sig)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
    # Index directly by name to preserve the exact row order
    sig_mat <- as.matrix(intensity[sig_ids_list, , drop = FALSE])
    return(sig_mat)
    }

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <-function(de_intensity, num_colors, palette) {
    colors <- brewer.pal(num_colors, palette)
    heatmap(
        de_intensity,
        col     = colors,
        margins = c(8, 4),
        cexCol  = 0.6,
        labRow  = NA
    )
}
