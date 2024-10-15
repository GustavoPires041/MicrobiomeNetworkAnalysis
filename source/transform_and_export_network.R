#' Transform and Export Network
#'
#' This function merges bacterial and fungal phyloseq objects, extracts necessary data,
#' calculates abundance, transforms taxonomic information, constructs a network, and
#' exports taxonomic and network data to specified files.
#'
#' @param phyloseq_bact Phyloseq object for bacteria
#' @param phyloseq_fung Phyloseq object for fungi
#' @param spiec_out SPiecEasi output
#' @param path_net Path to export the network file
#' @param path_tax Path to export the taxonomic file
#'
#' @return An igraph object representing the constructed network
#' @export
transform_and_export_network <- function(phyloseq_bact, phyloseq_fung, spiec_out, path_net, path_tax){
  library(SpiecEasi)
  library(igraph)
  
  # Merge phyloseq objects for bacteria and fungi
  merged_phyloseq <- merge_phyloseq(phyloseq_bact, phyloseq_fung)
  
  # Extract OTU table data
  otu.c <- t(otu_table(merged_phyloseq)@.Data)
  
  # Extract taxonomic information
  tax.c <- as.data.frame(tax_table(merged_phyloseq)@.Data)
  
  # Calculate abundance and transform to log2 percentages
  tax.c$Abundance <- taxa_sums(merged_phyloseq)
  
  # Replace spaces with underscores in the taxonomic column
  tax.c[, 6] <- gsub(" ", "_", tax.c[, 6])
  
  # Obtain the adjacency matrix from the SPiecEasi output
  n.c <- symBeta(getOptBeta(spiec_out))
  
  # Set column and row names for the adjacency matrix
  colnames(n.c) <- rownames(n.c) <- colnames(otu.c)
  
  # Create an undirected graph from the adjacency matrix
  net <- graph.adjacency(n.c, mode = 'undirected', add.rownames = TRUE, weighted = TRUE)
  
  # Write taxonomic information to a text file
  write.table(tax.c, file = path_tax, sep = "\t", quote = FALSE)
  
  # Write the network to a text file in ncol format
  write.graph(net, file = path_net, format = "ncol")
  
  # Return the constructed igraph
  return(net)
}
