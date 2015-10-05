#' Parse generic Level 3 TCGA data files
#' 
#' The TCGA Level 3 data files seem to all follow a similar
#' pattern with the first row being the sample names and the
#' second row having value descriptions ("Beta_value", 
#' "Expression", etc.). This does the first-level parsing
#' to give back the "data" as a data.frame and the first
#' row as column names. The second row is returned as the
#' second item in the list.
#' 
#' @note This function is not meant to be called by an
#' end-user, but it is exported for the purposes of extensibility
#' by a developer.
#' 
#' @param uri a url or filename of a TCGA data matrix
#' 
#' @return A list with elements \code{val}, a tbl_df and
#' \code{header} with the second row of the header as a 
#' character vector.
#' 
#' @export
#' 
#' @importFrom readr read_delim
readGenericTCGA = function(uri) {
  header = read_tsv(uri, n_max = 2, col_names = FALSE)
  vals    = read_tsv(uri, skip = 2, col_names = FALSE)
  colnames(vals) = as.character(header[1, ])
  return(list(vals = vals, header = as.character(header[2, ])))
}

#' Read Level 3 TCGA gene expression data
#' 
#' This function reads a file of Level 3 gene expression
#' data and returns a basic \code{SummarizedExperiment}
#' object with the data in \code{assay(object, 'exprs')}.
#' 
#' @param uri a url or filename of a TCGA data matrix
#' 
#' @import SummarizedExperiment
#' 
#' @export
#' 
#' @return A \code{SummarizedExperiment} object with
#' values in the \code{assay(object, 'exprs')},
#' sample names as colnames, and identifiers as 
#' rownames.
#' 
readGETCGA = function(uri) {
  tmp = readGenericTCGA(uri)
  se = SummarizedExperiment(assays = list(exprs = as.matrix(dplyr::select(tmp$vals,-1))))
  rownames(se) = make.unique(tmp$vals[[1]])
  return(se)
}

#' Read Level 3 TCGA methulation data
#' 
#' This function reads a file of Level 3 methylation data
#' data and returns a basic \code{SummarizedExperiment}
#' object with the data in \code{assay(object, 'exprs')}.
#' 
#' @param uri a url or filename of a TCGA data matrix
#' 
#' @import SummarizedExperiment
#' 
#' @export
#' 
#' @return A \code{SummarizedExperiment} object with
#' values in the \code{assay(object, 'beta')},
#' sample names as colnames, and identifiers as 
#' rownames.
#' 
readMethTCGA = function(uri) {
  tmp = readGenericTCGA(uri)
  # there are missing values in the mapping to the genome, so simply remove.
  tmp$vals = tmp$vals[!is.na(tmp$vals[,4]),]
  se  = SummarizedExperiment(assays = list(betas = as.matrix(tmp$vals[,tmp$header=='Beta_value'])),
                             rowRanges = GRanges(tmp$vals[[4]],
                                                 ranges = IRanges(start=tmp$vals[[5]],width=1),
                                                 Symbol = tmp$vals[,3]))
  rownames(se) = make.unique(tmp$vals[[1]])
  return(se)
}

#' Split TCGA IDs to pieces
#' 
#' @return A data frame with the portions of the TCGA ID
#' split into named columns.
#' 
#' @import stringi, dplyr
#' 
#' @export
#' 
splitTCGAIds = function(ids) {
  df = tbl_df(data.frame(do.call(rbind,stri_split_fixed(ids,pattern='-')),
                         stringsAsFactors = FALSE))
  colnames(df) = c('TCGA','TSS','Subject','SampleVial','PortionAnalyte','Plate','Center')
  return(df)
}


#' Write a matrix out as an Annotated Feature Matrix file
#' 
#' The Annotated Feature Matrix (.afm) format is described in 
#' detail at \url{https://code.google.com/p/rf-ace/wiki/Manual}.
#' Basically, though, features are in rows and samples/observations
#' are in columns. The rows are annotated with a data type that can 
#' include: 
#' \itemize {
#'   \item{B} for binary data
#'   \item{C} for categorical data (string)
#'   \item{N} for numeric data
#' }
#' 
#' @param mat an R matrix
#' @param fname the filename to which to write the resulting
#' afm file
#' @param append A logical, whether or not to append the matrix
#' to a file or overwrite the file. Note that columns should 
#' match the existing file if append=TRUE.
#' 
#' @details The column names are required as specified in
#' \url{https://code.google.com/p/rf-ace/wiki/Manual}. The row 
#' names are also required and will be written out with the 
#' appropriate prefix (B, N, C) based on the type of the matrix
#' to be written. 
#' 
#' @export
#' 
#' @examples
#' mat = matrix(rnorm(25),nc=5)
#' colnames(mat) = rownames(mat) = LETTERS[1:5]
#' tmpfile = tempfile()
#' write.afm(mat,tmpfile)
#' readLines(tmpfile)
write.afm = function(mat, fname, append=FALSE, ...) {
  prefix = "N"
  append = as.logical(append)
  if(is.character(mat)) prefix = "C"
  if(is.logical(mat)) prefix = "B"
  if(is.null(colnames(mat))) stop('column names are required')
  if(is.null(rownames(mat))) stop('row names are required')
  rnames = rownames(mat)
  rnames = paste(prefix, rnames, sep=":")
  write.table(data.frame(featureid=rnames,mat),
              file=fname,append = append, sep = "\t",
              col.names=!append, row.names=FALSE, quote=FALSE, ...)
}
