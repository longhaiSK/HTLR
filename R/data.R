#' Colon Tissues
#'
#' In this dataset, expression levels of 40 tumor and 22 normal colon tissues 
#' for 6500 human genes are measured using the Affymetrix technology. 
#' A selection of 2000 genes with highest minimal intensity across the samples
#' has been made by Alon et al. (1999). The data is preprocessed by carrying out 
#' a base 10 logarithmic transformation and standardizing each tissue sample to 
#' zero mean and unit variance across the genes.
#'
#' @docType data
#' 
#' @keywords datasets
#'
#' @format A list contains data matrix \code{X} and response vector \code{y}:
#' \describe{
#'   \item{X}{A matrix with 66 rows (observations) and 2000 columns (features).}
#'   \item{y}{A binary vector where 0 indicates normal colon tissues and 1 indicates tumor colon tissues.}
#' }
#' 
#' @usage data("colon")
#' 
#' @source \url{ftp://stat.ethz.ch/Manuscripts/dettling/bagboost.pdf}
#' 
#' @references Dettling Marcel, and Peter Bühlmann (2002). Supervised clustering of genes.
#' \emph{Genome biology}, 3(12), research0069-1.
#' 
"colon"

#' Pima Indians Diabetes
#'
#' This dataset is originally from the National Institute of Diabetes and Digestive and Kidney Diseases.
#' The objective of the dataset is to diagnostically predict whether or not a patient has diabetes, 
#' based on certain diagnostic measurements included in the dataset. Several constraints were placed 
#' on the selection of these instances from a larger database. In particular, all patients here are 
#' females at least 21 years old of Pima Indian heritage. Different from the UCI original version, 
#' the dataset has been preprocessed such that rows with missing values are removed, and features are scaled.
#'
#' @docType data
#' 
#' @keywords datasets
#'
#' @format A list contains data matrix \code{X} and response vector \code{y}:
#' \describe{
#'   \item{X}{A matrix with 392 rows (observations) and 8 columns (features).}
#'   \item{y}{A binary vector where 1 indicates diabetes patients and 0 for otherwise.}
#' }
#' 
#' @usage data("diabetes392")
#' 
#' @source \url{https://www.kaggle.com/uciml/pima-indians-diabetes-database}
#' 
#' @seealso \url{https://avehtari.github.io/modelselection/diabetes.html}
#' 
#' @references Smith, J.W., Everhart, J.E., Dickson, W.C., Knowler, W.C., & Johannes, R.S. (1988). 
#' Using the ADAP learning algorithm to forecast the onset of diabetes mellitus. 
#' \emph{In Proceedings of the Symposium on Computer Applications and Medical Care} (pp. 261--265). 
#' IEEE Computer Society Press.
#' 
"diabetes392"
