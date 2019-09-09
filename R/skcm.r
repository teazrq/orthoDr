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
#' <https://cancergenome.nih.gov/>
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
#' <http://bioinformatics.cing.ac.cy/MelGene/>
#' <https://cancergenome.nih.gov/>
#' @keywords skcm.melgene
"skcm.melgene"
