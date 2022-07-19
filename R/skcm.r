#' Skin Cutaneous Melanoma Data set
#' 
#' The clinical variables of the SKCM dataset. 
#' The original data was obtained from The Cancer Genome Atlas (TCGA).
#' 
#' @format 
#' Contains 469 subjects with 156 failures. Each row contains one subject,
#'  subject ID is indicated by row name. Variables include:
#'  
#' - `Time`
#' - `Censor`
#' - `Gender`
#' - `Age`
#' 
#' **Note:** `Age` has 8 missing values.
#' 
#' @references 
#' <https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga>
"skcm.clinical"

#' Genes associated with Melanoma given by the MelGene Database
#' 
#' The expression of top 20 genes of cutaneous melanoma literature based on the
#' MelGene Database.
#' @format
#' Each row contains one subject, subject ID is indicated by row name. 
#' Gene names in the columns. The columns are scaled.
#' 
#' @references 
#' Chatzinasiou, Foteini, Christina M. Lill, Katerina Kypreou, Irene Stefanaki, Vasiliki Nicolaou, George Spyrou, Evangelos Evangelou et al. "Comprehensive field synopsis and systematic meta-analyses of genetic association studies in cutaneous melanoma." Journal of the National Cancer Institute 103, no. 16 (2011): 1227-1235.
#' Emmanouil I. Athanasiadis, Kyriaki Antonopoulou, Foteini Chatzinasiou, Christina M. Lill, Marilena M. Bourdakou, Argiris Sakellariou, Katerina Kypreou, Irene Stefanaki, Evangelos Evangelou, John P.A. Ioannidis, Lars Bertram, Alexander J. Stratigos, George M. Spyrou, A Web-based database of genetic association studies in cutaneous melanoma enhanced with network-driven data exploration tools, Database, Volume 2014, 2014, bau101, https://doi.org/10.1093/database/bau101
#' <https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga>
#' @keywords skcm.melgene
"skcm.melgene"
