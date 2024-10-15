#' Identifica as taxas principais em um objeto phyloseq
#'
#' Esta função identifica as taxas principais com base em cortes de prevalência e abundância em um objeto phyloseq.
#'
#' @param pseq Objeto phyloseq contendo os dados de sequenciamento.
#' @param corte_prev O corte de prevalência para identificar as taxas principais (padrão: 0).
#' @param corte_abund O corte de abundância para identificar as taxas principais (padrão: 0.1).
#' @return Um data frame contendo informações das taxas principais identificadas.
#' @examples
#' \dontrun{
#' library(microbiome)
#' data("enterotype")
#' core_bacteria(enterotype, corte_prev = 0.1, corte_abund = 0.05)
#' }
#' @export
core = function(pseq, corte_prev, corte_abund){
  library(microbiome)
  pseq <- microbiome::transform(pseq, transform = "compositional")
  # Identifica as taxas principais com base nos cortes de prevalência e abundância
  core_taxa <- core_members(pseq, detection = corte_abund, prevalence = corte_prev)
  # Verifica se foram encontradas taxas principais
  if (length(core_taxa) == 0) {
    print("Nenhuma taxa principal encontrada com os cortes especificados.")
    return(NULL)
  }
  
  # Obtém a tabela de taxonomia do objeto phyloseq
  taxonomy <- tax_table(pseq)
  
  # Filtra as taxas principais na tabela de taxonomia
  core_taxa_info <- as.data.frame(taxonomy[rownames(taxonomy) %in% core_taxa, ])
  cat("Numero de membros restantes após filtro:",length(rownames(core_taxa_info)),"\n")
  # Retorna as informações das taxas principais
  return(core_taxa_info)
}
