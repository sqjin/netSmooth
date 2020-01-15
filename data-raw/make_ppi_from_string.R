require(STRINGdb)
require(igraph)
require(biomaRt)

getSTRINGdbForSpecies <- function(species=c('human','mouse')) {
  species <- match.arg(species)
  if(species=='human') return( STRINGdb$new(species=9606) )
  if(species=='mouse') return( STRINGdb$new(species=10090) )
  
}

getBiomartForspecies <- function(species=c('human','mouse')) {
  species <- match.arg(species)
  if(species=='human') return( useMart(host = 'grch37.ensembl.org',
                                       biomart='ENSEMBL_MART_ENSEMBL',
                                       dataset='hsapiens_gene_ensembl') )
  if(species=='mouse') return( useMart("ensembl",
                                       dataset="mmusculus_gene_ensembl") )
}

getPPIFromStringDB <- function(species) {
  string_db <- getSTRINGdbForSpecies(species)
  human_graph <- string_db$get_graph()
  
  # edge.scores <- E(human_graph)$combined_score
  # ninetyth.percentile <- quantile(edge.scores, 0.9)
  # thresh <- data.frame(name='90th percentile', val=ninetyth.percentile)
  # ninetyth.percentile <- 900
  #human_graph <- subgraph.edges(human_graph,E(human_graph)[combined_score >ninetyth.percentile])
  # Columns: protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
  # human_graph <- subgraph.edges(human_graph,E(human_graph)[database >= 900])
   human_graph <- subgraph.edges(human_graph,E(human_graph)[experiments >= 900]) # good choice
   #human_graph <- subgraph.edges(human_graph,E(human_graph)[combined_score >= 700])
  
  adj_matrix <- as_adjacency_matrix(human_graph)
  
  protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                        function(x) x[2])
  
  mart = getBiomartForspecies(species)
  
  # mart_results <- getBM(attributes = c("ensembl_gene_id","ensembl_peptide_id"),
  #                       filters = "ensembl_peptide_id", values = protein_ids,mart = mart)
  # mart_results <- getBM(attributes = c("external_gene_name","mgi_symbol","ensembl_gene_id","ensembl_peptide_id"),
  #                       filters = "ensembl_peptide_id", values = protein_ids,mart = mart)
  if (species=='human') {
    mart_results <- getBM(attributes = c("hgnc_symbol","ensembl_peptide_id"),
                          filters = "ensembl_peptide_id", values = protein_ids,mart = mart)
  } else if (species=='mouse') {
    mart_results <- getBM(attributes = c("mgi_symbol","ensembl_peptide_id"),
                          filters = "ensembl_peptide_id", values = protein_ids,mart = mart)
  }
  colnames(mart_results) <- c("gene_sysmbol","ensembl_peptide_id")
  ix <- match(protein_ids, mart_results$ensembl_peptide_id)
  ix <- ix[!is.na(ix)]
  ix <- intersect(ix, which(mart_results$gene_sysmbol != ""))
  
  
  newnames <- protein_ids
  newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
    mart_results[ix, 'gene_sysmbol']
  rownames(adj_matrix) <- newnames
  colnames(adj_matrix) <- newnames
  
  ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
  nullrows <- Matrix::rowSums(ppi)==0
  ppi <- ppi[!nullrows,!nullrows]
  return(ppi)
}

human.ppi <- getPPIFromStringDB('human')
mouse.ppi <- getPPIFromStringDB('mouse')

devtools::use_data(human.ppi, overwrite = TRUE)
devtools::use_data(mouse.ppi, overwrite = TRUE)
