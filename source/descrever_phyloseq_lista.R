#' Calcula estatísticas resumidas de abundância
#'
#' Esta função recebe um vetor de abundância e calcula várias estatísticas
#' resumidas, como número de elementos, mínimo, primeiro quartil, mediana,
#' terceiro quartil, máximo e média.
#'
#' @param abundance Vetor numérico contendo os valores de abundância
#' @return Um vetor com as estatísticas resumidas: número de elementos, mínimo,
#' primeiro quartil, mediana, terceiro quartil, máximo e média.
#'
#' @examples
#' abundance <- c(1, 2, 3, 4, 5)
#' calcular_estatisticas(abundance)
#'
calcular_estatisticas <- function(abundance) {
  # Calcula o número de elementos no vetor
  # e cria um vetor com as estatísticas desejadas
  stats <- c(
    taxa = length(abundance),
    minimo = min(abundance),
    primeiro_quartil = quantile(abundance, 0.25),
    mediana = median(abundance),
    terceiro_quartil = quantile(abundance, 0.75),
    maximo = max(abundance),
    media = mean(abundance)
  )
  
  return(stats)  # Retorna o vetor de estatísticas
}

#' Compara phyloseqs e calcula estatísticas resumidas de abundância
#'
#' Esta função recebe dois objetos phyloseq, calcula a abundância para cada um
#' deles e, em seguida, calcula estatísticas resumidas para compará-los.
#'
#' @param physeq1 Primeiro objeto phyloseq
#' @param nome1 Nome associado ao primeiro objeto phyloseq
#' @param physeq2 Segundo objeto phyloseq
#' @param nome2 Nome associado ao segundo objeto phyloseq
#' @return Um data frame contendo as estatísticas resumidas para os dois objetos phyloseq
#'
#' @examples
#' # physeq1 e physeq2 são objetos phyloseq
#' comparar_phyloseq(physeq1, "Amostra1", physeq2, "Amostra2")
#'
comparar_phyloseq <- function(physeq1, nome1, physeq2, nome2) {
  # Calcula a abundância para cada objeto phyloseq
  abundance1 <- taxa_sums(physeq1)
  abundance2 <- taxa_sums(physeq2)
  
  # Calcula as estatísticas resumidas para cada objeto phyloseq
  stats_physeq1 <- calcular_estatisticas(abundance1)
  stats_physeq2 <- calcular_estatisticas(abundance2)
  
  # Retorna um data frame com as estatísticas resumidas para os dois objetos phyloseq
  return(data.frame(nome1 = stats_physeq1, nome2 = stats_physeq2))
}

#' Descreve e compara os valores de abundância em uma lista de phyloseqs
#'
#' Esta função recebe uma lista de objetos phyloseq e compara os valores de
#' abundância entre eles, gerando um data frame com as estatísticas resumidas
#' para cada comparação.
#'
#' @param lista_phyloseqs Lista contendo objetos phyloseq
#' @return Um data frame com as estatísticas resumidas das comparações de abundância
#'
#' @examples
#' # lista_phyloseqs é uma lista de objetos phyloseq
#' descrever_phyloseqs_lista(lista_phyloseqs)
#'
descrever_phyloseqs_lista <- function(lista_phyloseqs) {
  num_phyloseqs <- length(lista_phyloseqs)
  
  # Cria uma matriz para armazenar os resultados das comparações
  resultados <- matrix(NA, nrow = 7, ncol = num_phyloseqs)
  rownames(resultados) <- c("Taxa","Min", "Q1", "Median", "Q3", "Max", "Mean")
  colnames(resultados) <- names(lista_phyloseqs)
  
  # Loop para comparar cada par de objetos phyloseq na lista
  for (i in 1:(num_phyloseqs - 1)) {
    for (j in (i + 1):num_phyloseqs) {
      # Compara os objetos phyloseq e obtém as estatísticas resumidas
      df_stats <- comparar_phyloseq(lista_phyloseqs[[i]], names(lista_phyloseqs)[i], lista_phyloseqs[[j]], names(lista_phyloseqs)[j])
      
      # Armazena as estatísticas na matriz de resultados
      resultados[, i] <- df_stats[, 1]
      resultados[, j] <- df_stats[, 2]
    }
  }
  
  # Retorna um data frame com as estatísticas resumidas das comparações
  return(as.data.frame(resultados))
}
