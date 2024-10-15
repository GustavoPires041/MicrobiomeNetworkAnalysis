#' Realiza a limpeza e organização de um objeto phyloseq de dados de microbioma
#'
#' Esta função executa uma série de etapas de limpeza em um objeto phyloseq contendo dados de microbioma.
#' Ela renomeia colunas, remove informações específicas, ajusta metadados e filtra os dados.
#'
#' @param fungi Objeto phyloseq contendo os dados de microbioma.
#' @return Um objeto phyloseq resultante após as operações de limpeza e organização.
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("Your_Phyloseq_Object")  # Substitua 'Your_Phyloseq_Object' com seus dados reais
#' clean_pseq(Your_Phyloseq_Object)
#' }
#' @export
clean_pseq <- function(fungi){
  library(phyloseq)  # Carrega a biblioteca phyloseq para trabalhar com dados de microbioma
  library(Biostrings)  # Carrega a biblioteca Biostrings para operações com strings
  
  # Renomeia as colunas da tabela de taxonomia para os níveis taxonômicos padrão
  colnames(tax_table(fungi)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Remove os prefixos específicos dos identificadores taxonômicos na tabela de taxonomia
  tax_table(fungi)[, colnames(tax_table(fungi))] <- gsub(tax_table(fungi)[, colnames(tax_table(fungi))], pattern = "D_[0-7]__", replacement = "")
  
  # Remove linhagens específicas (como Mitocôndrias e Cloroplastos) da tabela de taxonomia
  fungi <- subset_taxa(fungi, Family != "Mitochondria" & Order != "Chloroplast")
  
  # Ajusta coluna SampleType_Sites, removendo informações redundantes sobre a saúde
  fungi@sam_data$SampleType_Sites <- gsub("_Healthy", "", fungi@sam_data$SampleType_Sites)
  fungi@sam_data$SampleType_Sites <- gsub("_SD", "", fungi@sam_data$SampleType_Sites)
  
  # Cria coluna SampleReads com o número de reads nos metadados
  fungi@sam_data$SampleReads <- sample_sums(fungi)
  
  set.seed(1412)  # Define uma semente para a geração de números aleatórios
  
  # Calcula a cobertura dos dados de microbioma
  cov <- metagMisc::phyloseq_coverage(fungi)
  
  # Obtém os metadados do objeto phyloseq
  meta <- fungi@sam_data
  meta$SampleID <- rownames(meta)
  
  # Combina os dados de metadados com a cobertura calculada
  a <- merge(data.frame(meta), as.data.frame(cov), by = "SampleID")
  rownames(a) <- a$SampleID
  a$SampleID <- NULL
  a$SampleCoverage.x <- NULL
  a$SampleCoverage.y <- NULL
  
  # Atualiza os dados de metadados no objeto phyloseq
  fungi@sam_data <- sample_data(a)
  
  # Verifica se o primeiro elemento da tabela taxonômica é "Bacteria" e realiza ações específicas
  # Remove prefixos taxonômicos em potencial e filtra dados de acordo com condições
  for(i in c("k__", "p__", "c__", "o__", "f__", "g__", "s__")){
    tax_table(fungi)[, colnames(tax_table(fungi))] <- gsub(tax_table(fungi)[, colnames(tax_table(fungi))], pattern = i, replacement = "")
  }
  fungi <- subset_taxa(fungi, fungi@tax_table[, 1] != "No blast hit")
  fungi.f <- subset_samples(fungi, SampleReads > 2500 & SampleCoverage > 0.97)
  
  # Separa o objeto phyloseq com base em variáveis específicas
  split1 <- metagMisc::phyloseq_sep_variable(fungi.f, "SampleType_Sites")
  split2 <- metagMisc::phyloseq_sep_variable(fungi.f, "SampleType_Surface_Object")
  ps.split <- split1[c("Floor", "Wall", "Scalp", "Hand")]
  ps.split$Object <- split2$Object
  
  return(ps.split)  # Retorna o objeto phyloseq resultante após as operações
}
