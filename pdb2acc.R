# pdb2acc.R

################ PDB_ID to ACC #############################################
# This function converts a PDB identifier into an UniProt accession number.

pdb2acc <- function(pdbID, email){
  
  require(httr)
  
  uniprot_url <- "http://www.uniprot.org/uploadlists/"
  
  params <- list(from = 'PDB_ID', 
                 to = 'ACC',
                 format = 'tab',
                 query = pdbID)
  
  my_headers <- add_headers('User-Agent' = paste('R', email))
  
  request <- GET(uniprot_url, query = params, my_headers)
  
  ans <- content(request, 'text', encoding = "ISO-8859-1")
  
  ## ------- Ancillary Función ----------
  extract_acc <- function(acc_id, ans){
    acc_list <- NULL
    len_acc <- nchar(acc_id)
    pos <- gregexpr(acc_id, ans)[[1]]
    for (i in 1:length(pos)){
      pos_tab <- pos[i] + gregexpr('\t', substr(ans, pos[i], nchar(ans)))[[1]] - 1
      id <- substr(ans, pos_tab[2]+1, pos_tab[3]-1)
      acc_list <- c(acc_list, id)
    }
    acc_list
  }
  ## ------------------------------------
  
  return(extract_acc(pdbID, ans))
}
#######################################################################################