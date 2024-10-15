#' Realiza a limpeza e organização de um objeto phyloseq de dados de microbioma
#'
#' Esta função executa uma série de etapas de limpeza em um objeto phyloseq contendo dados de microbioma.
#' Ela renomeia colunas, remove informações específicas, ajusta metadados e filtra os dados.
#'
#' @param bacteria Objeto phyloseq contendo os dados de microbioma.
#' @return Um objeto phyloseq resultante após as operações de limpeza e organização.
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("Your_Phyloseq_Object")  # Substitua 'Your_Phyloseq_Object' com seus dados reais
#' clean_bacteria(Your_Phyloseq_Object)
#' }
#' @importFrom phyloseq sample_sums subset_samples tax_table sample_data
#' @importFrom Biostrings gsub
#' @export
clean_bacteria <- function(bacteria){
  library(phyloseq)  # Carrega a biblioteca phyloseq para trabalhar com dados de microbioma
  library(Biostrings)  # Carrega a biblioteca Biostrings para operações com strings
  
  # Renomeia as colunas da tabela de taxonomia para os níveis taxonômicos padrão
  colnames(tax_table(bacteria)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Remove os prefixos específicos dos identificadores taxonômicos na tabela de taxonomia
  tax_table(bacteria)[, colnames(tax_table(bacteria))] <- gsub(tax_table(bacteria)[, colnames(tax_table(bacteria))], pattern = "D_[0-7]__", replacement = "")
  
  # Remove linhagens específicas (como Mitocôndrias e Cloroplastos) da tabela de taxonomia
  bacteria <- subset_taxa(bacteria, Family != "Mitochondria" & Order != "Chloroplast")
  
  # Ajusta coluna SampleType_Sites, removendo informações redundantes sobre a saúde
  bacteria@sam_data$SampleType_Sites <- gsub("_Healthy", "", bacteria@sam_data$SampleType_Sites)
  bacteria@sam_data$SampleType_Sites <- gsub("_SD", "", bacteria@sam_data$SampleType_Sites)
  
  # Cria coluna SampleReads com o número de reads nos metadados
  bacteria@sam_data$SampleReads <- sample_sums(bacteria)
  
  set.seed(1412)  # Define uma semente para a geração de números aleatórios
  
  # Calcula a cobertura dos dados de microbioma
  cov <- metagMisc::phyloseq_coverage(bacteria)
  
  # Obtém os metadados do objeto phyloseq
  meta <- bacteria@sam_data
  meta$SampleID <- rownames(meta)
  
  # Combina os dados de metadados com a cobertura calculada
  a <- merge(data.frame(meta), as.data.frame(cov), by = "SampleID")
  rownames(a) <- a$SampleID
  a$SampleID <- NULL
  a$SampleCoverage.x <- NULL
  a$SampleCoverage.y <- NULL
  
  # Atualiza os dados de metadados no objeto phyloseq
  bacteria@sam_data <- sample_data(a)
  
  # Filtra os dados com base na contagem de reads e cobertura
  bacteria.f <- subset_samples(bacteria, SampleReads > 1000 & SampleCoverage > 0.97)
  
  # Separa o objeto phyloseq com base em variáveis específicas
  split1 <- metagMisc::phyloseq_sep_variable(bacteria.f, "SampleType_Sites")
  split2 <- metagMisc::phyloseq_sep_variable(bacteria.f, "SampleType_Surface_Object")
  ps.split <- split1[c("Floor", "Wall", "Scalp", "Hand")]
  ps.split$Object <- split2$Object
  
  return(ps.split)  # Retorna o objeto phyloseq resultante após as operações
}

