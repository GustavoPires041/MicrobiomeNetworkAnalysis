##Selecao do core
core_fungi = function(fungi, corte_prev = 0, corte_abund){
  fungi = fungi
  corte_prev = corte_prev
  corte_abund = corte_abund
  ## Função para criar um core de taxa em um objeto fungioseq
  
  ## Transforma o objeto fungioseq para um formato composicional
  fungi_compositional <- microbiome::transform(fungi, "compositional") 
  
  ## Identifica as taxas principais com base nos cortes de abundância e prevalência
  core_taxa <- core_members(fungi_compositional, detection = corte_abund, prevalence = corte_prev)
  print(core_taxa)
  ## Obtém a tabela de taxonomia das taxas identificadas como principais
  taxonomy <- as.data.frame(tax_table(fungi_compositional))
  core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core_taxa)
  
  ## Exibe informações do objeto fungioseq original
  print(fungi)
  
  ## Cria um novo objeto fungioseq contendo apenas as taxas principais identificadas
  fungi_core = subset_taxa(fungi, taxa_names(fungi) %in% core_taxa) 
  
  ## Retorna o novo objeto fungioseq com as taxas principais
  return(fungi_core)
}