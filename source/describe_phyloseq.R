#' Imprime informações descritivas de um objeto phyloseq e plota a distribuição da abundância total dos taxa.
#'
#' Esta função recebe um objeto phyloseq, extrai informações como número de taxas, número de amostras
#' e estatísticas de resumo da abundância total dos taxa. Além disso, plota um gráfico de densidade
#' representando a distribuição da abundância total dos taxa no conjunto de dados.
#'
#' @param physeq Um objeto phyloseq.
#' @param nome O nome atribuído ao objeto phyloseq.
#'
#' @return Esta função não retorna um objeto específico, apenas imprime informações e plota um gráfico.
#'
#' @import ggplot2
#' @importFrom scales label_number
#'
#' @examples
#' describe_phyloseq(my_phyloseq_object, "MeuObjetoPhyloseq")
#'
#' @export
describe_phyloseq <- function(physeq, nome) {
  library(ggplot2)  # Carrega a biblioteca ggplot2 para criação de gráficos
  
  # Obtém os nomes das taxas e das amostras
  taxa <- taxa_names(physeq)
  samples <- sample_names(physeq)
  
  # Calcula a abundância total dos taxa por amostra
  abundance <- sample_sums(physeq)
  
  # Imprime o nome do elemento, número de taxa e número de amostras
  cat("\033[1m\nNome do Elemento:\033[0m", nome, "\n")  # Imprime o nome em negrito
  cat("Número de Taxa:", length(taxa), "\n")  # Imprime o número de taxa
  cat("Número de Amostras:", length(samples), "\n")  # Imprime o número de amostras
  
  # Calcula as estatísticas de resumo da abundância
  summary_stats <- summary(abundance)
  
  # Imprime as estatísticas de resumo da abundância
  cat("\033[1m\nEstatísticas de Abundância:\033[0m\n")  # Título em negrito
  cat("Mínimo:", summary_stats[1], "\n")  # Imprime o valor mínimo
  cat("Primeiro Quartil (Q1):", summary_stats[2], "\n")  # Imprime o primeiro quartil
  cat("Mediana (Q2):", summary_stats[3], "\n")  # Imprime a mediana
  cat("Terceiro Quartil (Q3):", summary_stats[4], "\n")  # Imprime o terceiro quartil
  cat("Máximo:", summary_stats[5], "\n")  # Imprime o valor máximo
  cat("Média:", mean(abundance), "\n\n")  # Imprime a média
  
  # Cria um dataframe para o gráfico de densidade
  data <- data.frame(Samples = names(abundance), Abundance = abundance)
  
  # Cria o gráfico de distribuição da abundância total dos taxa (density plot)
  p <- ggplot(data, aes(x = Abundance)) +
    geom_density(fill = "skyblue", color = "black") +  # Adiciona a camada de densidade
    labs(title = paste("Distribuição da Abundância Total dos Taxa -", nome), x = "Abundância Total", y = "Frequência") +  # Define os títulos dos eixos
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, accuracy = 1))  # Formatação do eixo y
  
  # Mostra o gráfico
  print(p)
}
