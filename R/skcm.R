#' @title Clinical Data for Cutaneous Melanoma
#' @name skcm.clinical
#' @description Clinical variables for the skin cutaneous melanoma (SKCM)
#'   data from The Cancer Genome Atlas (TCGA).
#' @usage skcm.clinical
#' @docType data
#' @format Contains 469 subjects with 156 events. Each row represents one
#'   subject, identified by its row name. Variables are \code{Time},
#'   \code{Censor}, \code{Gender}, and \code{Age}. \code{Age} has 8 missing values.
#' @references \url{https://www.cancer.gov/ccg/research/genome-sequencing/tcga}
#' @keywords skcm.clinical
NULL

#' @title Melanoma Gene-Expression Data
#' @name skcm.melgene
#' @description Expression measurements for 20 genes selected from the
#'   cutaneous melanoma literature using the MelGene Database.
#' @usage skcm.melgene
#' @docType data
#' @format Each row represents one subject, identified by its row name.
#'   Columns are scaled gene-expression measurements, with gene names as
#'   column names.
#' @references Chatzinasiou, Foteini, Christina M. Lill, Katerina Kypreou, Irene Stefanaki, Vasiliki Nicolaou, George Spyrou, Evangelos Evangelou et al. "Comprehensive field synopsis and systematic meta-analyses of genetic association studies in cutaneous melanoma." Journal of the National Cancer Institute 103, no. 16 (2011): 1227-1235.
#' @references \url{https://www.cancer.gov/ccg/research/genome-sequencing/tcga}
#' @keywords skcm.melgene
NULL

