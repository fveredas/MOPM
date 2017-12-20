# dif.R

############  DIFFERENCES BETWEEN PDB AND UNIPROT ENUMERATIONS  #######################
# Sometimes the PDB and UniProt proteins are enumerated in different ways.
# This function found the difference 'dif' in the following way: Uniprot = PDB + dif

dif <- function(pdbID, chain, email){
  
  # UniProt
  ACC <- pdb2acc(pdbID, email)[1]
  download.file(paste("http://www.uniprot.org/uniprot/", ACC, ".fasta", sep=""), destfile="./temp.fas")
  uniprot <- read.fasta("./temp.fas")[[2]]
  system("rm temp.fas")
  uniprot <- toupper(paste(uniprot, collapse=""))
  
  # PDB
  mypdb <- bio3d::read.pdb(pdbID)
  pdb.seq <- aa321(mypdb$atom$resid[which(mypdb$atom$elety == "CA" & mypdb$atom$type == "ATOM" &
                                            mypdb$atom$chain == chain)])
  
  pos <- mypdb$atom[which(mypdb$atom$elety == "CA"  & mypdb$atom$type == "ATOM" &
                            mypdb$atom$chain == chain),]$resno
  names(pdb.seq) <- pos
  probe <- pdb.seq[10:20]
  probe <- paste(probe, collapse="") # from pos 10 to 20 in the pdb
  
  # Hybridization using the pdb probe
  d <- (regexpr(probe, uniprot)[[1]]) - pos[10]
  return(d)
}
#######################################################################################