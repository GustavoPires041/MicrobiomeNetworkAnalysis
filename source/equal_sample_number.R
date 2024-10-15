#' Equalize sample numbers within subgroups before splitting phyloseq objects
#'
#' This function equalizes the number of samples within subgroups ("Healthy" and "SD") 
#' for two phyloseq objects before splitting based on the specified category.
#'
#' @param phyloseq_bac A phyloseq object for the first dataset.
#' @param phyloseq_fun A phyloseq object for the second dataset.
#' @param category The variable used for splitting the phyloseq objects (default: "HouseStatus").
#'
#' @return A list containing equalized subgroups for both phyloseq objects: bac_healthy, bac_sd, 
#'         fun_healthy, fun_sd.
#'
#' @details This function aims to equalize the number of samples within each subgroup 
#'          ("Healthy" and "SD") for both phyloseq objects before performing a split 
#'          based on the specified category. It uses the prune_samples function to 
#'          ensure an equal number of samples in each subgroup before splitting.
#'
#' @examples
#' # Usage example
#' result <- equal_sample_number(phyloseq_data1, phyloseq_data2, category = "HouseStatus")
#'
#' @export
equal_sample_number <- function(phyloseq_bac, phyloseq_fun, category = "HouseStatus") {
  # Get sample names for both phyloseq objects
  sample_names_bac <- sample_names(phyloseq_bac)
  sample_names_fun <- sample_names(phyloseq_fun)
  
  # Find common sample names
  common_samples <- intersect(sample_names_bac, sample_names_fun)
  
  # Subset phyloseq objects to keep only common samples
  phyloseq_bac <- prune_samples(common_samples, phyloseq_bac)
  phyloseq_fun <- prune_samples(common_samples, phyloseq_fun)
  
  # Split into subgroups
  phyloseq_bac_split <- metagMisc::phyloseq_sep_variable(phyloseq_bac, category)
  phyloseq_fun_split <- metagMisc::phyloseq_sep_variable(phyloseq_fun, category)
  
  # Equalize number of samples within each subgroup for both phyloseq objects
  for (subgroup in c("Healthy", "SD")) {
    # Get number of samples in each subgroup for both phyloseq objects
    n_samples_bac <- nsamples(phyloseq_bac_split[[subgroup]])
    n_samples_fun <- nsamples(phyloseq_fun_split[[subgroup]])
    
    # Find common samples within the subgroup
    common_samples_subgroup <- intersect(sample_names(phyloseq_bac_split[[subgroup]]), 
                                         sample_names(phyloseq_fun_split[[subgroup]]))
    
    # Equalize number of samples
    if (n_samples_bac > n_samples_fun) {
      phyloseq_bac_split[[subgroup]] <- prune_samples(common_samples_subgroup[1:n_samples_fun], 
                                                      phyloseq_bac_split[[subgroup]])
    } else if (n_samples_bac < n_samples_fun) {
      phyloseq_fun_split[[subgroup]] <- prune_samples(common_samples_subgroup[1:n_samples_bac], 
                                                      phyloseq_fun_split[[subgroup]])
    }
  }
  
  # Return the equalized subgroups for both phyloseq objects
  return(c(bac_healthy = phyloseq_bac_split$Healthy,
           bac_sd = phyloseq_bac_split$SD,
           fun_healthy = phyloseq_fun_split$Healthy,
           fun_sd = phyloseq_fun_split$SD))
}
