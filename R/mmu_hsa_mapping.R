#' Mouse to Human Gene Mapping
#'
#' A data frame containing gene id mappings between mouse (Mus musculus) genes andhuman (Homo sapiens) genes.
#'
#' @format A data frame with 6 columns
#' \describe{
#'   \item{mmu_symbol}{mouse gene names}
#'   \item{mmu_ensembl}{mouse gene ensembl ID}
#'   \item{mmu_entrez}{mouse gene entrez ID}
#'   \item{hsa_symbol}{human gene names}
#'   \item{hsa_ensembl}{human gene ensembl ID}
#'   \item{hsa_entrez}{human gene entrez id}
#' }
#' @source bioMart
#' @name mmu_hsa_mapping
"mmu_hsa_mapping"