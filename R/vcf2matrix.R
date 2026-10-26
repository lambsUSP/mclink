#' Read a vcf file and returns a list with the chromossomes, the positions and the genotypes
#'
#' @param path path to vcf file
#' @param chr_num number of the chromosome of interest (optional)
#' @param ref genomic reference to use
#' @param limit memory limit to use to read the file
#' @return list with the chromossomes, the positions of the SNPs and a sparse matrix with the genotypes of the vcf file
#' @importFrom vcfR read.vcfR getCHROM getPOS
#' @importFrom Matrix Matrix
#' @export
vcf2matrix <- function(path, chr_num = 0, ref = NULL, limit = NULL) {

  vcf <- if(!is.null(limit)) read.vcfR(path, limit = limit) else read.vcfR(path)

  if (chr_num != 0) vcf <- get_chromosome(vcf, chr_num)

  chrs <- getCHROM(vcf)
  
  if(!is.null(ref)){
    keep <- startsWith(chrs, ref)
    vcf <- vcf[keep, ]
  }

  chrs <- chrs |> factor() |> as.integer()

  pos <- getPOS(vcf)

  geno <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
  geno[is.na(geno)] <- -1
  dimension <- dim(geno)
  geno <- Matrix(data = geno, sparse = TRUE)

  list(geno, chrs, pos)
}


#' Get the chromosome of interest.
#' Assume that the vcf file is sorted by chromosome.
#'
#' @param vcf vcf file
#' @param chr_num number of the chromoseme of interest
#' @return vcf object with just the chromosome chosen
#' @importFrom vcfR extract.gt getCHROM
get_chromosome <- function(vcf, chr_num) {
  chr_unique <- getCHROM(vcf) |> unique() 
  chromosome <- chr_unique[chr_num]

  if (is.na(chromosome)) paste("Chromosome ", chr_num, " out of the range. Chromosomes numbers avaliable are: ", paste(seq_along(chr_unique), collapse = ", ")) |> stop()
  idx <- vcf@fix[, "CHROM"] == chromosome
  vcf <- vcf[idx, ]

  vcf
}
