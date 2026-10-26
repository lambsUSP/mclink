#' Analyze genetic distances by chromosome
#'
#' Splits the SNP matrix by chromosome and runs \code{simulate_panmixia()} independently for each chromosome and for the complete genome.
#'
#' @param data         List returned by \code{vcf2matrix()}, containing the genotype matrix, chromosome indices, and SNP positions.
#' @param chr          Chromosome to be analyzed
#' @param iterations   Number of panmictic simulations. Default \code{1000}.
#' @param sample_size  Number of SNPs sampled per iteration. Default \code{1000}.
#' @param min_distance Minimum distance between sampled SNPs (bp). Default \code{1000}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{genome}{Results for the complete genome.}
#'     \item{chromosomes}{Named list with results per chromosome.}
#'   }
#' @export


genetic_distribution <- function(data, 
                             chr = NULL,
                             iterations = 1000,
                             sample_size = 1000,
                             min_distance = 1000) {
  
  geno <- data[[1]]
  chromosomes <- data[[2]]
  pos  <- data[[3]]
  
  chr_unique <- sort(unique(chromosomes))
  
  result_chr <- vector("list", length(chr_unique))
  names(result_chr) <- paste0("Chromosome ", sprintf("%02d", chr_unique))

  if(is.null(chr)){
    result_genome <- simulate_panmixia(data,
                                        iterations   = iterations,
                                        sample_size  = sample_size,
                                        min_distance = min_distance)
    result_chr = NULL 
  }else if(chr == "by_chr"){
    for (it_chr in chr_unique) {
      cat("Processing chromosome", it_chr, "\n")
      idx <- which(chromosomes == it_chr)

      sub_data <- list(
                       geno[idx, ],
                       chromosomes[idx],
                       pos[idx]
     )

      result_chr[[it_chr]] <- simulate_panmixia(sub_data, 
                                           iterations = iterations, 
                                           sample_size = sample_size,
                                           min_distance = min_distance)

    }

    cat("Processing whole genome\n")
    result_genome <- simulate_panmixia(data,
                                       iterations   = iterations,
                                       sample_size  = sample_size,
                                       min_distance = min_distance)
  }else{
      if (!(chr %in% chr_unique)) 
         paste("The chromosome ", chr," is not valid, the chromosomes numbers available are: ", paste(chr_unique, collapse = ", ")) |> stop()
      
      cat("Processing chromosome", chr, "\n")
      idx <- which(chromosomes == chr)
      sub_data <- list(
                       geno[idx, ],
                       chr[idx],
                       pos[idx]
      )
      result_chr[[chr]] <- simulate_panmixia(sub_data, 
                                           iterations = iterations, 
                                           sample_size = sample_size,
                                           min_distance = min_distance)
      result_genome = NULL
  }

  distributions <- list(genome      = result_genome,
       chromosomes = result_chr,
       sample_size = sample_size) 

  class(distributions) <- "mclink_genetic_distribution"

  distributions
}


#' Plot genetic distances got by genetic_distribution function 
#'
#' @param result object returned by genetic_distribution()
#' @param chr chromosome to be ploted
#' @param rd boolean to decide if should plot the genetic distance or the rd plot
#' @param ... extended parameters
#'
#' @method plot mclink_genetic_distribution
#' @export
plot.mclink_genetic_distribution <- function(result, chr = NULL, rd = FALSE, ...){
  if(!is.null(chr)){
    if (rd)
      plot_association_index(result, chr = chr, ...)
    else
      plot_genetic_distance(result, chr = chr, ...)
  }

  else if(!is.null(result$genome)){
    if (rd)
      plot_association_index(result, ...)
    else 
      plot_genetic_distance(result$genome, ...)
  }else{
    for (i in seq_along(result$chromosomes)) {
      if(!is.null(result$chromosomes[[i]])){
        if (rd)
          plot_association_index(result, chr = i, ...)
        else
          plot_genetic_distance(result, chr = i, ...)

      } 
    }
  }
}
