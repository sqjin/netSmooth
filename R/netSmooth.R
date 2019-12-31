setGeneric(
    name = "netSmooth",
    def = function(x, ...) {
        standardGeneric("netSmooth")
    }
)

#' Perform network smoothing of gene expression or other omics data
#' @param x    matrix or single cell Experiment
#' @param adjMatrix    adjacency matrix of gene network to use
#' @param alpha    numeric in [0,1]
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @param filepath      String: Path to location where hdf5 output file is supposed to be saved.
#'                      Will be ignored when regular matrices or SummarizedExperiment are
#'                      used as input.
#' @return network-smoothed gene expression matrix or SummarizedExperiment
#'         object
#' @examples
#' x <- matrix(rnbinom(12000, size=1, prob = .1), ncol=60)
#' rownames(x) <- paste0('gene', seq_len(dim(x)[1]))
#'
#' adj_matrix <- matrix(as.numeric(rnorm(200*200)>.8), ncol=200)
#' rownames(adj_matrix) <- colnames(adj_matrix) <- paste0('gene', seq_len(dim(x)[1]))
#' x.smoothed <- netSmooth(x, adj_matrix, alpha=0.5)
#' @export
#' @rdname netSmooth
#' @inheritParams netSmooth,matrix-method
#' @aliases netSmooth
setMethod("netSmooth",
    signature(x='matrix'),
    function(x, adjMatrix, alpha=0.5,
        normalizeAdjMatrix=c('rows','columns')
        ) {
        normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)

        stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
        stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)))
        if(sum(Matrix::rowSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
        if(sum(Matrix::colSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")

        if(is.numeric(alpha)) {
            if(alpha<0 | alpha > 1) {
                stop('alpha must be between 0 and 1')
            }
            x.smoothed <- smoothAndRecombine(x, adjMatrix, alpha,
                                        normalizeAdjMatrix=normalizeAdjMatrix)
        } else stop("unsupported alpha value: ", class(alpha))
        return(x.smoothed)
    }
)
#' @param x    single cell Experiment
#' @param adjMatrix    adjacency matrix of gene network to use
#' @param alpha    numeric in [0,1]
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
#' @rdname netSmooth
#' @export
setMethod("netSmooth",
          signature(x='SingleCellExperiment'),
          function(x, adjMatrix, alpha=0.5,
                   normalizeAdjMatrix=c('rows','columns')) {
            matrixdata <- assay(x)
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            ret <- netSmooth(matrixdata, adjMatrix, alpha,
                             normalizeAdjMatrix)
            return(SingleCellExperiment(assays = list(counts = ret)))
          })

#' @rdname netSmooth
#' @export
setMethod("netSmooth",
          signature(x='Matrix'),
          function(x, adjMatrix, alpha=0.5,
                   normalizeAdjMatrix=c('rows','columns')) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)

            stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
            stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)))
            if(sum(Matrix::rowSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
            if(sum(Matrix::colSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")

            if(is.numeric(alpha)) {
              if(alpha<0 | alpha > 1) {
                stop('alpha must be between 0 and 1')
              }
              x.smoothed <- smoothAndRecombine(x, adjMatrix, alpha,
                                               normalizeAdjMatrix=normalizeAdjMatrix)
            } else stop("unsupported alpha value: ", class(alpha))
            return(x.smoothed)
          }
)

#' @rdname netSmooth
#' @export
setMethod("netSmooth",
          signature(x='DelayedMatrix'),

          function(x, adjMatrix, alpha=0.5,
                   normalizeAdjMatrix=c('rows','columns'),
                   filepath = NULL)
          {

            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)

            stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
            stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)))
            if(sum(Matrix::rowSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")
            if(sum(Matrix::colSums(adjMatrix)==0)>0) stop("PPI cannot have zero rows/columns")


            if(is.numeric(alpha)) {
              if(alpha<0 | alpha > 1) {
                stop('alpha must be between 0 and 1')
              }

              x.smoothed <- smoothAndRecombine(x, adjMatrix, alpha,
                                               normalizeAdjMatrix=normalizeAdjMatrix,
                                               filepath=filepath)
            } else stop("unsupported alpha value: ", class(alpha))


            return(x.smoothed)
          })

