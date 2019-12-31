
#' Mouse Protein-Protein interaction graph
#'
#' An adjacency matrix of the high confidence known experimental interactions 
#' between mouse proteins on STRINGdb.
#'
#' See the script in `system.file(package="netSmooth", "data-raw", "make_ppi_from_string.R")`
#' for full details of how this object was made.
#'
#' @format A square matrix where A_{ij}=1 if gene i interacts with gene j
#' @source \url{http://www.string-db.org/}
"mouse.ppi"

#' Human Protein-Protein interaction graph
#'
#' An adjacency matrix of the high confidence known experimental interactions
#' between human proteins on STRINGdb.
#'
#' See the script in `system.file(package="netSmooth", "data-raw", "make_ppi_from_string.R")`
#' for full details of how this object was made.
#'
#' @format A square matrix where A_{ij}=1 if gene i interacts with gene j
#' @source \url{http://www.string-db.org/}
"human.ppi"
