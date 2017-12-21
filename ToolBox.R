# ToolBox.R

# Collection of R functions whose purpose is to assist us when computing features 
# related to physico-chemical and biological properties of methionyl residues.

# ------------------------------------------------------------------------------ #
# INTER-RESIDUE DISTANCES .................................. residueDist
# SAMPLE SEQ FROM FASTA .................................... sample.seq
# SEARCH AMPHIPATHIC ....................................... searchAmphipathic
# MUTATE METHIONINE TO X ................................... m2x
# HYDROPHOBIC MOMENT ....................................... hydrophobic.moment
# PARALOGOUS FILTER ........................................ pFilter
# RELABELING TIPS FROM A TREE .............................. rlabeling.tree
# RELABELING SPECIES FROM A MSA ............................ rlabeling.aln
# alnID .................................................... alnID
# PDB_ID to ACC ........................................... pdb2acc
# DIFFERENCES BETWEEN PDB AND UNIPROT ENUMERATIONS ......... dif
# B FACTOR FOR THE DELTA SULFUR ATOM FROM METIONINES ....... factorB
# SOLVENT ACCESSIBLE SURFACE AREA .......................... accDSSP
# ATOM DEPTH ANALYSIS ...................................... dpx
# EIGHT SPECIES ENTROPY ANALYSIS ........................... eightShannon
# SHANNON  ENTROPY ANALYSIS ................................ shannon
# S-AROMATIC MOTIF ANALYSIS ................................ SaromaticDistnace
# PHOSPHOSITES ANALYSIS .................................... phosphosite
# ------------------------------------------------------------------------------ #

library(bio3d)
library(httr)

############ INTER-RESIDUE DISTANCES  ####################################
# Copyright (c) 2017 Juan Carlos Aledo (MIT License)
# ----------------------------------------------------------------------------------
# --------------- Brief summary of what this sofware does --------------------------
# ----------------------------------------------------------------------------------
# The current program takes as input  a PDB identifier, the position
# of two residue from that protein and their corresponding chains. In 
# addition there are two logical arguments. One is 'backbone', 
# when TRUE it means that we include those atoms belonging to the main
# chain (CA, N, O and C) beside all the side chain atoms. The other is
# 'H.atoms', when TRUE we include all the hydrogen atoms in the 
# as long as the PDB provides their coordinates. This function returns
# a list of three elements, where each of these elements is, in turn, a
# list of three elments. 
#
# Example:
# myresult <- residueDist("1Q8K", 51, "A", 55, "A", backbone=T, H.atoms=T)
#
# myresult[[1]]
# [[1]]
# [1] 2.964445      --This is the minimum iteratomic distance, in angstroms, between:
# [[2]]
# [1] "SER-A-HB2"   -- atom HB2 from Serine at 51 in chain A and
# [[3]]
# [1] "ILE-A-HD13   --atom HD13 from Isoleucine at 55 in chain A.
#
# myresult[[2]]
# [[1]]
# [1] 9.278745      --This is the maximum iteratomic distance, in angstroms, between:
# [[2]]
# [1] "SER-A-HA"   -- atom HA from Serine at 51 in chain A and
# [[3]]
# [1] "ILE-A-O"   --atom O from Isoleucine at 55 in chain A.
#
# myresult[[3]]
# [[1]]
# [1] 5.961756      --This is the mean iteratomic distance, in angstroms, between:
# [[2]]
# [1] "SER-A"   --  Serine at 51 in chain A and
# [[3]]
# [1] "ILE-A"   --  Isoleucine at 55 in chain A.
# ------------------------------------------------------------------------------------

residueDist <- function(pdbID, A, chainA, B, chainB, backbone, H.atoms){
  
  #---------------- Matricial Function to Compute Distances
  dist.binomial <- function(a, b, squared=TRUE, ...){ 
    an <- apply(a, 1, function(x) crossprod(x, x))
    bn <- apply(b, 1, function(x) crossprod(x, x))
    
    m <- length(an)
    n <- length(bn)
    
    an_bn <- matrix(rep(an, n), nrow=m) + matrix(rep(bn, m), nrow=m, byrow=T)
    d2 <- an_bn - 2 * tcrossprod(a, b)
    
    if(!squared){
      d2 <- sqrt(d2)
    }
    attr(d2, "squared") <- squared
    return(d2)
  }
  
  # ------------ Data Preparation and Computation of Distances
  library(bio3d)
  mypdb <- bio3d::read.pdb(pdbID)
  atomsA <- mypdb$atom[which(mypdb$atom$resno == A & mypdb$atom$chain == chainA),]
  atA <- matrix(t(c(atomsA$x, atomsA$y, atomsA$z)), ncol=3)
  rownames(atA) <- atomsA$elety
  colnames(atA) <- c("x", "y", "z")
  atomsB <- mypdb$atom[which(mypdb$atom$resno == B & mypdb$atom$chain == chainB),]
  atB <- matrix(t(c(atomsB$x, atomsB$y, atomsB$z)), ncol=3)
  rownames(atB) <- atomsB$elety
  colnames(atB) <- c("x", "y", "z")
  
  bdist <- dist.binomial(atA, atB, squared = FALSE)
  
  #--------------- Main/Side Chain 
  if (backbone==FALSE){
    # Atoms belonging to the backbone are removed  
    bdist <- bdist[which(! rownames(bdist) %in% c("CA", "N", "C", "O")), 
                   which(! colnames(bdist) %in% c("CA", "N", "C", "O"))]
  } 
  
  #-------------- Hydrogen Atoms
  # Sometimes the pdb contains hydrogen atoms' coordinates
  if (H.atoms == FALSE){
    filas <- rownames(bdist)[which(substr(rownames(bdist), 1,1) != "H")]
    columnas <- colnames(bdist)[which(substr(colnames(bdist), 1,1) != "H")]
    bdist <- bdist[which(rownames(bdist) %in%  filas), 
                   which(colnames(bdist) %in% columnas)]
  }
  #--------------- Selecting Min, Max and Mean Distances
  min.d <- min(bdist)
  fila.min <- which(bdist==min(bdist), arr.ind=TRUE)[1]
  columna.min <- which(bdist==min(bdist), arr.ind=TRUE)[2]
  rA <- paste(atomsA$resid[1], "-", chainA, "-", rownames(bdist)[fila.min], sep="")
  rB <- paste(atomsB$resid[1], "-", chainB, "-", colnames(bdist)[columna.min], sep="")
  
  mylist <- list(list(min.d, rA, rB))
  
  max.d <- max(bdist)
  fila.max <- which(bdist==max(bdist), arr.ind=TRUE)[1]
  columna.max <- which(bdist==max(bdist), arr.ind=TRUE)[2]
  rA <- paste(atomsA$resid[1], "-", chainA, "-", rownames(bdist)[fila.max], sep="")
  rB <- paste(atomsB$resid[1], "-", chainB, "-", colnames(bdist)[columna.max], sep="")
  
  mylist[[2]] <- list(max.d, rA, rB)
  
  rA <- paste(atomsA$resid[1], "-", chainA, sep="")
  rB <- paste(atomsB$resid[1], "-", chainB, sep="")
  
  mylist[[3]] <- list(mean(bdist), rA, rB)
  
  return(mylist)
}
################################################################################

##################### SAMPLE SEQ FROM FASTA FILE ###############################
# This function returns the protein sequence of the requested species 
# (provided as an argument, 'species') extracted from the indicated 
# fasta file (also provided as an argument, 'fasFile').

# Example: 
# sample.seq(species = "Homo_sapiens", fasFile = "./CaM_selection.fas")

# Alternatively, the function can also manage fasta files provided 
# for eggNOG (where the species are identified by numbers).
# To this end, the argument 'eggNOT' must be turn to TRUE and the 
# user must supply an appropiated 'protein_taxon' argument (i.e. "CaM_Metazoa").

# Example:
# sample.seq(species = "Homo sapiens", eggNOG = T, protein_taxon = "CaM_Metazoa")

# In this last case, the files 'protein_taxon_list.txt' and 'protein_taxon.fas' 
# must be found in the working directory.

sample.seq <- function(species, eggNOG=F, protein_taxon, fasFile){

  if (eggNOG==T){
    mytable <- read.csv(file=paste("./", protein_taxon, "_list.txt", sep=""), header= F, sep="\t")
    allseq <- bio3d::read.fasta(file= paste("./", protein_taxon, ".fas", sep=""))
    
    lista <- strsplit(allseq$id, "\\.")
    lista <- lapply(lista, function(x) x[1])
    id <- mytable$V4[which(mytable$V3 == species)][1] # from latin name to id 
    t <- allseq$ali[which(lista == id), which(allseq$ali[which(lista==id),] != "-")]
    t <- paste(t, collapse="")
  }
  else {
    allseq <- bio3d::read.fasta(file=paste("./", fasFile, sep=""))
    sp1 <- gsub(" ", "_", species)
    sp2 <- gsub("_", " ", species)
    target <- c(species, sp1, sp2)
    index <- which(allseq$id %in% target)
    t <- allseq$ali[index, which(allseq$ali[index,] != "-")]
    t <- paste(t, collapse="")
  }
  
  return(t)
}
########################################################################

#############   SEARCH AMPHIPATHIC  ####################################
# This function takes two arguments:
# 'polypeptide' is a named polypeptide sequence character vector
# (the names correspond to the residue position in the protein),
# and 'n' is the number of residues of the searched helix (window's 
# size). The function returns a matrix where each row is the start
# (in the polypeptide sequence) of the considered helix and each
# column is the hydrophobic moment at the indicated (column name)
# degree.

searchAmphipathic <- function(polypeptide, n){
  m <- length(polypeptide)
  H.matrix <- matrix(data=NA, nrow=(m-n+1), ncol=9)
  colnames(H.matrix)<-c(90, 101.25, 112.5, 123.75, 135, 146.25,
                        157.5, 168.75, 180)
  
  for (i in 1:(m-n+1)){
    h <- window(polypeptide, i, (i+n-1))
    H.matrix[i,] <- round(hydrophobic.moment(h),2)
  }
  return(H.matrix)
}
#################################################################################

############# MUTATE METHIONINE TO X ############################################
# This function takes the following arguments:
# (1) 'target' a named character vector containing an amino acid sequence
# (the names indicate the positions).
# (2) 'pos' a numeric vector containing the positions of the met to be mutated.
# (3) 'x' one letter code for the amino acid to be introduced instead Met.
# If no methionine is found at some of the indicated positions, the program
# will stop and complain. 
# if the argument 'pos' is absent, then all the methionines present in the 
# sequence will be changed.
# The function returns the mutated sequence as a character vector

m2x <- function(target, pos, x){
  if (!missing(pos)){
    # checking all the residues to be changed are methionines
    met <- target[which(names(target) %in% pos)]
    p <- which(met != "M")
    if (length(p) != 0){
      print(helice[names(p)])
      stop("\n IS/ARE NOT METHIONINE(S)")
    }
    else{
      mutant <- target
      mutant[names(mutant) %in% pos] <- x
    }
    return(mutant)
  }
  if (missing(pos)){
    mutant <- gsub("M", x, target)
    return(mutant)
  }
}
####################################################################################################

############# HYDROPHOBIC MOMENT  #################################################################
# This function takes an argument named 'segment' that is a vector containing the sequence of 
# amino acids that supposedly form a helix. An optional second argument is 'index' which can 
# takes two values, either 'ei' or 'bm' (Eisenberg or Black & Mould hydrophobicity indexes,
# respectively), the first one is taken by default if none is specified. The function returns a
# named numeric vector containing the hydrophobic moment a diferent angles (the angles in degrees
# are provided as names)

hydrophobic.moment <- function(segment, index="ei"){
  library(bio3d)
  Eisenberg <- aa.index[68][[1]]$I # vector con los valores de hidrofobicidad para cada Aa.
  BlackMould <- aa.index[528][[1]]$I # vector con los valores de hidrofobicidad para cada Aa.
  
  delta <- seq(from = pi/2, to = pi, by = pi/16) # conjunto de ángulos en radianes
  # Cada uno de estos ángulos son posibles valores para una estructura regular helicoidal.
  # Nótese que en el caso de hélices alfa, delta es 2pi/3.6 = 1.7279 rad (100 grados).
  # 'segment', el argumento de la función es un vector conteniendo los residuos de la hélice que se desea analizar.
  suma_sen <- 0
  suma_cos <- 0
  
  if (index %in% c("ei", "bm") == FALSE){
    cat(paste("  index = \"", index, "\" not found.\n", sep=""))
    cat("  Defaulting to index = \" ei \" (Eisenberg index) \n\n")
    index <- "ei"
  }
  
  if (index == "ei"){ 
    for (i in 1:length(segment)){
      Hn <- Eisenberg[which(names(Eisenberg) == segment[i])]
      suma_sen <- suma_sen + ( Hn*sin(delta*i) )
      suma_cos <- suma_cos + ( Hn*cos(delta*i) )
    }
    mu <- sqrt( (suma_sen)^2 + (suma_cos)^2 )
  }
  else if (index == "bm" ){
    for (i in 1:length(segment)){
      Hn <- BlackMould[which(names(BlackMould) == segment[i])]
      suma_sen <- suma_sen + ( Hn*sin(delta*i) )
      suma_cos <- suma_cos + ( Hn*cos(delta*i) )
    }
    mu <- sqrt( (suma_sen)^2 + (suma_cos)^2 )
  }
  else { stop("Error while computing the hydrophobic moment")}
  
  delta_g <- (180/pi)*seq(from = pi/2, to = pi, by = pi/16) # Pasamos de rad a grados 
  names(mu) <- delta_g
  return(mu)
}  
####################################################################################################

############# PARALOGOUS FILTER  ##################################################################
# This function takes as arguments a raw MSA that may contain paralogous as well as orthologus sequences, 
# and the uniprot ID of the reference sequence. The function returns the filtered alignment containing 
# only orthologous. When a species has different paralogous, the function takes the closer one to the 
# reference sequence. 

pFilter <- function(msa, uniprotID){
  require(bio3d)
  
  download.file(paste("http://uniprot.org/uniprot/", uniprotID, ".fasta", sep=""), destfile="./temp.fas")
  needle <- bio3d::read.fasta("./temp.fas")
  needle <- paste(needle$ali, collapse="") # Reference protein sequence
  system("rm temp.fas") # Unix
  # system("cmd.exe", input = "del temp.fas") # Windows
  
  seqs <- c()   
  raw_aln <- bio3d::read.fasta(msa) 
  for (i in 1:nrow(raw_aln$ali)){
    seqs <- c(seqs, paste(raw_aln$ali[i,], collapse="")) 
  }
  names(seqs) <- raw_aln$id
  
  homo_id_index <- which(regexpr("9606", raw_aln$id) != -1) # index of the reference sequence in raw_aln
  haystack <- seqs[homo_id_index] # Reference protein sequence containing gaps after alignment
  min_adist <- which.min(adist(needle, haystack, partial=T)) # Cual es la menor...
  
  # Empezamos generando un data frame donde vamos a ir metiendo, para cada especie, el idenfificador y la secuencia más cercana a la primera en la tabla, que será la de referencia en humanos:
  species_orthologous <- data.frame(species=names(haystack[min_adist]), sequence=haystack[min_adist])
  species_ids <- c("9606")
  
  
  # Vamos recorriendo los ortólogos
  for(i in 1:length(seqs)){
    id <- unlist(strsplit(names(seqs)[i], split='[.]'))[1]
    if (!(id %in% species_ids)){
      species <- which(regexpr(id, names(seqs)) == 1)
      haystack <- seqs[species] # Las secuencias de proteínas de esa especie alineadas
      ids <- names(haystack) # Los identificadores
      min_adist <- which.min(adist(needle, haystack, partial=T)) # Cual es la menor...
      species_orthologous <- rbind(species_orthologous, data.frame(species=ids[min_adist], sequence=haystack[min_adist]))
      species_ids <- c(species_ids, id)
    }
  }
  species_orthologous$sequence <- as.character(species_orthologous$sequence)
  
  filtered_aln <- c() # alignment once the paralogous have been removed
  for (i in 1:nrow(species_orthologous)){
    other <- strsplit(species_orthologous$sequence[i], split="")[[1]]
    filtered_aln <- seqbind(filtered_aln, other)
  }
  filtered_aln$id <- as.character(species_orthologous$species)
  return(filtered_aln)
}
#########################################################################################

############# RELABELING TIPS FROM A TREE ###############################################
# Renaming the tip label in the MSA from eggNOG. 
# For instance, the label for human (the human sequence)
# is "9606.ENSP00000256383", but I want it to apper just as "Homo_sapiens"

rlabeling.tree <- function(table, tree){
  for (i in 1:nrow(table)){
    target <- paste(table$V4[i], ".", table$V1[i], sep="")
    ## Changing the labels in the tree
    for (j in 1:length(tree$tip.label)){
      if (tree$tip.label[j] == target){
        print(as.character(table$V3[i]))
        tree$tip.label[j] <- as.character(table$V3[i])
      }
    }
  }
  return(tree)
}
#########################################################################################

##### RELABELING SPECIES FROM A MSA #####################################################
# Renaming the  the species label in the MSA from eggNOG. 
# For instance, the label for human (the human sequence)
# is "9606.ENSP00000256383", but I want it to apper just as "Homo_sapiens"

rlabeling.aln <- function(table, aln){
  for (i in 1:nrow(table)){
    target <- paste(table$V4[i], ".", table$V1[i], sep="")
    ## Changing the labels in the aln
    for (j in 1:length(aln$id)){
      if (as.character(aln$id[j]) == target){
        aln$id[j] <- as.character(table$V3[i])
      }
    }
  }
  return(aln)
}
########################################################################################

#############  alnID  ###################################################################
# This function takes as argument 'ids' a character vector containing the identifiers of
# the proteins whose sequences we want to align. An optional argument, 'sps' can be passed 
# to this function to indicates the origen (ie. species) of the sequences.

alnID <- function(ids, sps){
  
  # ------- Function to download the sequences ---------------------------------------- #
  descarga <- function(target){
    download.file(url=paste("http://www.uniprot.org/uniprot/", target, ".fasta", sep=""), 
                  destfile=paste("./", target, ".fas", sep=""))
    prot <-bio3d::read.fasta(paste("./", target, ".fas", sep=""))
    
    system(paste("rm ", target, ".fas", sep="")) # Don't keep the fasta files # Unix
    # system("cmd.exe", input = paste("del ", target, ".fas", sep="")) # Don't keep the fasta files # Windows
    
    return(prot)
  }
  # ---------------------------------------------------------------------------------- #
  
  # ------ Preparing a matrix containing the sequences 
  P <- lapply(ids, descarga)
  P <- lapply(P, function(x) x$ali)
  
  M <- P[[1]] 
  for (i in 2:length(P)){
    M <- bio3d::seqbind(M, P[[i]], blank="-")
  }
  # ---------------------------------------------------------------------------------- #
  
  
  # --------- The alinement is carried out by MUSCLE
  # z <- list(...)
  # if(!is.null(z$sps)) {
  #   aln <- seqaln(M, id=sps, exefile = "/Users/JCA/Dropbox/MySOFTWARE/MUSCLE/muscle" )
  # }
  # else {
  #   aln <- seqaln(M, exefile = "/Users/JCA/Dropbox/MySOFTWARE/MUSCLE/muscle")
  # }
  if (missing(sps)){
    # aln <- seqaln(M, exefile = "/Users/JCA/Dropbox/MySOFTWARE/MUSCLE/muscle")
    aln <- seqaln(M, exefile = "muscle")
  }
  else {
    # aln <- seqaln(M, id=sps, exefile = "/Users/JCA/Dropbox/MySOFTWARE/MUSCLE/muscle")
    aln <- seqaln(M, id=sps, exefile = "muscle")
  }
  # ---------------------------------------------------------------------------------- #
  
  return(aln)
}
#######################################################################################

################ PDB_ID to ACC ########################################################
# The following two functions work together to convert a PDB identifier into an 
# UniProt accession number.

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
  extract_acc(pdbID, ans)
}
#######################################################################################

############  DIFFERENCES BETWEEN PDB AND UNIPROT ENUMERATIONS  #######################
# Sometimes the PDB and UniProt proteins are enumerated in different ways.
# This function found the difference 'dif' such as Uniprot = PDB + dif

dif <- function(pdbID, chain, email){
  
  # UniProt
  ACC <- pdb2acc(pdbID, email)[1]
  download.file(paste("http://www.uniprot.org/uniprot/", ACC, ".fasta", sep=""), destfile="./temp.fas")
  uniprot <- read.fasta("./temp.fas")[[2]]
  system("rm temp.fas") # Unix
  # system("cmd.exe", input = "del temp.fas") # Windows
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

################   B FACTOR FOR THE DELTA SULFUR ATOM FROM METIONINES  ################
factorB <- function(pdbID){
  mypdb <- bio3d::read.pdb(pdbID)
  Bfactor <- mypdb$atom[which(mypdb$ato$elety == "SD"),]
  return(Bfactor)
}
#######################################################################################


#################   SOLVENT ACCESSIBLE SURFACE AREA    ###############################
# Nota:
# SASA es una de las características computadas para desarrollar los modelos predictivos
# descritos en BMC Bioinformatics. Para computar dicha característica se usó una versión
# distinta a la que se recoge más abajo y que se encuentra en el script "ToolBox.R" que
# se ha incluido en el paquete de ficheros proporcionado como material suplementario.

# A partir de un fichero PDB analizamos qué residues están expuestos en superficie y 
# cuáles se encuentran enterrados en el interior de la proteína. Esta función (accDSSP) 
# toma como argumentos el identificador PDB de la estructura de una proteína 'id' y la
# difrencia en enumeración entre la enumeración en el fichero pdb y la
# enumeración de la estructura primaria en un fichero fasta; cuando no exista discrepancia
# entre ambas enumeraciones, consignar dif = 0. El program devuelve tanto la supercie en
# angstroms expueta al solvente (sasa) como la fracción de accesibilidad (acc), la cual 
# se calcula dividiendo el valor de sasa calculado por la superficie que dicho aminoácido
# tendría en un tripéptido tripéptido GXG con el esqueleto polipeptídico en una 
# conformación extendida y la conformación de la cadena lateral la más frecuentemente 
# observada en las proteínas (J. Mol. Biol. 196: 641-656).

# Update 29 Sept 2017 Al correr esta función con id = "1CLL" dio error
# problema con los iones CA, se introducen modificaciones en tres líneas
# que están marcadas como MODIFICADO. Además, la opción met = F, nos da
# la posibilidad de obtener la accesibilidad de todos los residuos y no sólo de Met.

# Update 23 Oct 2017. Para que devuelva la accesibilidad de cada aminoácido
# tanto cuando se analiza el complejo como cuando se analiza la subunidad de 
# forma aislada (Esta última siempre ha de ser mayor)

# Update 30 Oct 2017. Para que el fichero pdb se guarde/se lea en/desde el paso
# actual en el que se esté trabajando. También para que las funciones apeladas
# que puedan estar enmascaradas sean del paquete deseado.

accDSSP <- function(id, dif, met = T, keepPDB = F){
  library(bio3d)
  pdbsplit(get.pdb(id)) # MODIFICADO 23-Oct-2017 
  # complex <- read.pdb(paste("/Users/JCA/Dropbox/Programacion/R/R_ToolBox/", id, ".pdb", sep=""))
  complex <- bio3d::read.pdb(paste( id, ".pdb", sep=""))  # MODIFICADO 30 OCT
  
  pdb.seq <- aa321(complex$atom$resid[which(complex$atom$elety == "CA")])
  pdb.seq <- pdb.seq[which(pdb.seq != "X")]
  
  cx <- dssp(complex)
  
  acc <- data.frame(Aa = pdb.seq)
  acc$ResidNr <- complex$atom$resno[which(complex$atom$elety == "CA")][which(pdb.seq != "X")] + dif ## MODIFICADO
  acc$Chain <- complex$atom$chain[which(complex$atom$elety == "CA")][which(pdb.seq != "X")] ## MODIFICADO
  
  acc$sasa <- cx$acc
  acc$acc <- -999
  acc$acc[which(acc$Aa == "A")] <- round(acc$sasa[which(acc$Aa == "A")]/113, 3)
  acc$acc[which(acc$Aa == "R")] <- round(acc$sasa[which(acc$Aa == "R")]/241, 3)
  acc$acc[which(acc$Aa == "N")] <- round(acc$sasa[which(acc$Aa == "N")]/158, 3)
  acc$acc[which(acc$Aa == "D")] <- round(acc$sasa[which(acc$Aa == "D")]/151, 3)
  acc$acc[which(acc$Aa == "C")] <- round(acc$sasa[which(acc$Aa == "C")]/140, 3)
  acc$acc[which(acc$Aa == "Q")] <- round(acc$sasa[which(acc$Aa == "Q")]/189, 3)
  acc$acc[which(acc$Aa == "E")] <- round(acc$sasa[which(acc$Aa == "E")]/183, 3)
  acc$acc[which(acc$Aa == "G")] <- round(acc$sasa[which(acc$Aa == "G")]/85, 3)
  acc$acc[which(acc$Aa == "H")] <- round(acc$sasa[which(acc$Aa == "H")]/194, 3)
  acc$acc[which(acc$Aa == "I")] <- round(acc$sasa[which(acc$Aa == "I")]/182, 3)
  acc$acc[which(acc$Aa == "L")] <- round(acc$sasa[which(acc$Aa == "L")]/180, 3)
  acc$acc[which(acc$Aa == "K")] <- round(acc$sasa[which(acc$Aa == "K")]/211, 3)
  acc$acc[which(acc$Aa == "M")] <- round(acc$sasa[which(acc$Aa == "M")]/204, 3)
  acc$acc[which(acc$Aa == "F")] <- round(acc$sasa[which(acc$Aa == "F")]/218, 3)
  acc$acc[which(acc$Aa == "P")] <- round(acc$sasa[which(acc$Aa == "P")]/143, 3)
  acc$acc[which(acc$Aa == "S")] <- round(acc$sasa[which(acc$Aa == "S")]/122, 3)
  acc$acc[which(acc$Aa == "T")] <- round(acc$sasa[which(acc$Aa == "T")]/146, 3)
  acc$acc[which(acc$Aa == "W")] <- round(acc$sasa[which(acc$Aa == "W")]/259, 3)
  acc$acc[which(acc$Aa == "Y")] <- round(acc$sasa[which(acc$Aa == "Y")]/229, 3)
  acc$acc[which(acc$Aa == "V")] <- round(acc$sasa[which(acc$Aa == "V")]/160, 3)
  
  names(acc)[4] <- "sasa.complex"
  names(acc)[5] <- "acc.complex"
  
  ## ---- Isolated subunits 
  cadenas <- unique(complex$atom$chain)
  sacc <- sub.seq <- c()
  for (i in 1:length(cadenas)){
    current.chain <- cadenas[i]
    path <- paste("./split_chain/", id, "_", current.chain, ".pdb", sep="")
    single <- bio3d::read.pdb(path) # MODIFICADO 30 Oct, añadido el bio3d::
    sx <- dssp(single)
    sacc <- c(sacc, sx$acc)
    sub.seq <- c(sub.seq, aa321(single$atom$resid[which(single$atom$elety == "CA")]))
  }
  
  if (sum(sub.seq == acc$Aa)/nrow(acc) == 1) {
    acc$sasa.subunit <- sacc
    acc$acc.subunit <- -999
    acc$acc.subunit[which(acc$Aa == "A")] <- round(acc$sasa.subunit[which(acc$Aa == "A")]/113, 3)
    acc$acc.subunit[which(acc$Aa == "R")] <- round(acc$sasa.subunit[which(acc$Aa == "R")]/241, 3)
    acc$acc.subunit[which(acc$Aa == "N")] <- round(acc$sasa.subunit[which(acc$Aa == "N")]/158, 3)
    acc$acc.subunit[which(acc$Aa == "D")] <- round(acc$sasa.subunit[which(acc$Aa == "D")]/151, 3)
    acc$acc.subunit[which(acc$Aa == "C")] <- round(acc$sasa.subunit[which(acc$Aa == "C")]/140, 3)
    acc$acc.subunit[which(acc$Aa == "Q")] <- round(acc$sasa.subunit[which(acc$Aa == "Q")]/189, 3)
    acc$acc.subunit[which(acc$Aa == "E")] <- round(acc$sasa.subunit[which(acc$Aa == "E")]/183, 3)
    acc$acc.subunit[which(acc$Aa == "G")] <- round(acc$sasa.subunit[which(acc$Aa == "G")]/85, 3)
    acc$acc.subunit[which(acc$Aa == "H")] <- round(acc$sasa.subunit[which(acc$Aa == "H")]/194, 3)
    acc$acc.subunit[which(acc$Aa == "I")] <- round(acc$sasa.subunit[which(acc$Aa == "I")]/182, 3)
    acc$acc.subunit[which(acc$Aa == "L")] <- round(acc$sasa.subunit[which(acc$Aa == "L")]/180, 3)
    acc$acc.subunit[which(acc$Aa == "K")] <- round(acc$sasa.subunit[which(acc$Aa == "K")]/211, 3)
    acc$acc.subunit[which(acc$Aa == "M")] <- round(acc$sasa.subunit[which(acc$Aa == "M")]/204, 3)
    acc$acc.subunit[which(acc$Aa == "F")] <- round(acc$sasa.subunit[which(acc$Aa == "F")]/218, 3)
    acc$acc.subunit[which(acc$Aa == "P")] <- round(acc$sasa.subunit[which(acc$Aa == "P")]/143, 3)
    acc$acc.subunit[which(acc$Aa == "S")] <- round(acc$sasa.subunit[which(acc$Aa == "S")]/122, 3)
    acc$acc.subunit[which(acc$Aa == "T")] <- round(acc$sasa.subunit[which(acc$Aa == "T")]/146, 3)
    acc$acc.subunit[which(acc$Aa == "W")] <- round(acc$sasa.subunit[which(acc$Aa == "W")]/259, 3)
    acc$acc.subunit[which(acc$Aa == "Y")] <- round(acc$sasa.subunit[which(acc$Aa == "Y")]/229, 3)
    acc$acc.subunit[which(acc$Aa == "V")] <- round(acc$sasa.subunit[which(acc$Aa == "V")]/160, 3)
    
  } else { stop("Check number residues")}
  
  acc$DeltaSASA <- NA # sasa.subunit - sasa.complex
  for(i in 1:nrow(acc)){
    acc$DeltaSASA[i] <- acc$sasa.subunit[i] - acc$sasa.complex[i]
  }
  
  
  if (keepPDB == F) {
    system(paste("rm ", id, ".pdb", sep="")) # Unix
    system("rm -rf split_chain") # Unix
    # system("cmd.exe", input = paste("del ", id, ".pdb", sep="")) # Windows
    # system("cmd.exe", input = "rd /s /q split_chain") # Windows
  }
  
  
  if (met == F){
    return(acc)
  } else if (met == T){
    acc <- acc[which(acc$Aa == "M"),]
  } else {
    stop("The met argument must be a boolean")
  }
}
#########################################################################################

###################     ATOM DEPTH ANALYSIS      ########################################
# This function computes the sulfur delta atom's depth for each methionine residue

dpx <- function(pdbID, email){
  
  library(bio3d)
  get.pdb(pdbID)
  mypdb <- bio3d::read.pdb(pdbID)
  ath <- mypdb$atom # All atoms including heteroatoms
  at <- ath[which(ath$type == "ATOM"),] # Removing heteroatoms
  
  ############### Making use of GetArea
  require(RCurl)
  
  # Parameters
  url <- "http://curie.utmb.edu/cgi-bin/getarea.cgi"
  my_pdbs <- pdbID
  water <- "1.4"
  gradient <- "n"
  name <- "test"
  email <- email
  Method <- "4"
  
  result <- postForm(url, 
                     "water" = water, 
                     "gradient" = gradient, 
                     "name" = name, 
                     "email" = email, 
                     "Method" = Method, 
                     "PDBfile" = fileUpload(paste(my_pdbs, ".pdb", sep="")))
  
  writeLines(gsub("[</pre></td>|<td><pre>]","", result), con=paste(my_pdbs, ".txt", sep=''))
  res <- read.table(paste("./", pdbID, ".txt",sep=""), header=TRUE, as.is=TRUE,row.names = NULL, fill = TRUE, skip=16)
  res <- res[1:(nrow(res)-15),]
  
  at$acc <- NA
  for (i in 1:nrow(at)){
    for (j in i:nrow(res)){
      if(at$elety[i] == res$ATOM[j] & at$resno[i] == res$RESIDUE[j]){
        at$acc[i] <- res$AREAENERGY[j]
        break
      }
    }
  }
  at$acc <- as.numeric(at$acc)
  
  buried <- at[which(at$acc == 0),]
  exposed <- at[which(at$acc != 0),]
  buried$dpx <- NA # By defect. It will be changed later on when required.
  
  dist.binomial <- function(a, b, squared=TRUE, ...){ # Function to compute distances
      an <- apply(a, 1, function(x) crossprod(x, x))
      bn <- apply(b, 1, function(x) crossprod(x, x))

      m <- length(an)
      n <- length(bn)

      an_bn <- matrix(rep(an, n), nrow=m) + matrix(rep(bn, m), nrow=m, byrow=T)
      d2 <- an_bn - 2 * tcrossprod(a, b)

      if(!squared){
          d2 <- sqrt(d2)
      }
      attr(d2, "squared") <- squared
      return(d2)
  }
 
  points_buried <- matrix(t(c(buried$x, buried$y, buried$z)), ncol = 3)
  points_exposed <- matrix(t(c(exposed$x, exposed$y, exposed$z)), ncol = 3)
  bdist <- dist.binomial(points_buried, points_exposed, squared = FALSE)
  buried$dpx <- apply(bdist, 1, min)
  
  SDexp <- exposed[which(exposed$elety == "SD"),] # SD:  Sulfur Delta
  if (nrow(SDexp) != 0) { # Otherwise, it gives problem when there is not exposed SD atoms
    SDexp$dpx <- 0
  }
  SDbur <- buried[which(buried$elety == "SD"),]
  SD <- rbind(SDexp, SDbur)
  SD <- SD[with(SD, order(chain, resno)),]
  SD <- SD[,c(6,7,18)] # chain, resno and dpx
  
  system(paste("rm", paste("./", pdbID, ".txt",sep=""), sep=" ")) # Unix
  # system("cmd.exe", input = paste("del", paste(pdbID, ".txt",sep=""), sep=" ")) # Windows
  
  return(SD)
}
#########################################################################################

######### EIGHT SPECIES ENTROPY ANALYSIS ################################################
# Besides the target sequence for the protein under study (target) for which we need the 
# Uniprot accession (ACC), the orthologous proteins from  Pan troglodytes,
# Gorilla gorilla, Rattus norvegicus, Bos taurus, Gallus gallus, 
# Xenopus tropicalis and Danio rerio are aligned. This alignment is used to compute the 
# Shannon entropy with the help of the function 'shannon' from the current toolbox.

eightShannon <- function(target, ACC){
  
  # Loading the databases that will be required
  load("./btaurus_orth.rda")
  load("./drerio_orth.rda")
  load("./ggallus_orth.rda")
  load("./ggorilla_orth.rda")
  load("./ptroglodytes_orth.rda")
  load("./rnorvegicus_orth.rda")
  load("./xtropicalis_orth.rda")
  
  if (ACC == "P0DMV8"){ ACC <- "P08107"} # Exception for this ID
  chimp <- ptroglodytes_orth$SEQUENCE[which(ptroglodytes_orth$ACC_ID == ACC)]
  gorilla <- ggorilla_orth$SEQUENCE[which(ggorilla_orth$ACC_ID == ACC)]
  rat <- rnorvegicus_orth$SEQUENCE[which(rnorvegicus_orth$ACC_ID == ACC)]
  cow <- btaurus_orth$SEQUENCE[which(btaurus_orth$ACC_ID == ACC)]
  chicken <- ggallus_orth$SEQUENCE[which(ggallus_orth$ACC_ID == ACC)]
  frog <- xtropicalis_orth$SEQUENCE[which(xtropicalis_orth$ACC_ID == ACC)]
  fish <- drerio_orth$SEQUENCE[which(drerio_orth$ACC_ID == ACC)]
  if (ACC == "P08107") { ACC <- "P0DMV8"} # Restored.
  
  # Check that there exist orthologous sequences for all the 8 species
  if (length(target) == 0 | length(chimp) == 0 | length(gorilla) == 0 | length(rat) == 0 |
      length(cow) == 0 | length(chicken) == 0 | length(frog) == 0 | length(fish) == 0){ 
    output.Shannon <- NA
  } else if (target != "" & chimp != "" & gorilla != "" & rat != "" & cow != "" & 
      chicken != "" & frog != "" & fish != ""){
    
    seqs <- seqbind(target, chimp, gorilla, rat, cow, chicken, frog, fish, blank="-")
    vert <- seqaln(seqs, id= c("target", "chimp", "gorilla","rat","cow","chicken","frog", "fish"), exefile="muscle")
    system("rm aln.fa") # Unix
    # system("cmd.exe", input = "del aln.fa") # Windows
    
    # In the alignment, all the columns containing a gap in the target sequence are removed.
    # Since we are interested in those columns where a Met appears in the target sequence, this procedure
    # allows keeping the target sequence numeration without any detriment for our goals:
    vert$ali <- vert$ali[, which(vert$ali[1,] != "-")]
    
    output.Shannon <- shannon(vert)  # Now the Shannon entropy is computed.
    
  } else { # Cuando el alineamiento es incompleto (falta la secuencia de alguna especie).
    output.Shannon <- NA
  }
  return(output.Shannon)
}
#########################################################################################

################   SHANNON  ENTROPY ANALYSIS   ##########################################
shannon <- function(aln){ # La función toma un alineamiento como argumento
  myDF <- as.data.frame(aln$ali)
  # myDF<- as.data.frame(myDF[,which(myDF[1,]== aa)]) 
  # Es preciso forzar a dataframe cuando sólo hay una Met conservada en todas las especies.
  tablas <- lapply(myDF, table)
  tablas <- lapply(tablas, as.data.frame)
  
  # Por cada aa en la secuencia de referencia tenemos una $V. Ejemplo de la estructura de datos:
  #   $V1
  #   Var1 Freq
  #   1    -    2
  #   2    L    1
  #   3    M    5
  
  H <- function(df) { # for each position we compute the Shannon entropy
    -sum( ( df[,2]/(sum( df[,2] ) ) )*log(df[,2]/(sum(df[,2])), base=21) )
  }
  
  entropy <- sapply(tablas, H) # Vector with the Shannon entropies for each position
  entropy <- sapply(entropy, function(x) round(x, digits=3))
  Mean.entropy <- mean(entropy)
  SD.entropy <- sd(entropy)
  
  aa.pos <- which(myDF[1,] == "M") # Posiciones en las que aparece Met
  aa.H <- entropy[aa.pos] # Entropía en tales posiciones
  
  # To compute the frequency of Met in each relevant column of the alignment:
  temp <- aln$ali[,which(aln$ali[1,] == "M")] # Sub-alignment for Met position in the human protein
  aa.fM <- vector(mode = "numeric")
  for (i in 1:ncol(temp)){
    aa.fM[i] <- length(temp[,i][which(temp[,i] == "M")])
  }
  
  aa.fM <- aa.fM/nrow(temp)
  res <- list(aa.pos, aa.H, Mean.entropy, SD.entropy, aa.fM)
  return(res)
}
#################################################################################################


#####################    S-AROMATIC MOTIF ANALYSIS    ###########################################
# This function takes three arguments: (1) the pdb id or the full path to the pdb file. 
# (2) If we used the pdb id, the second argument should be online = TRUE (by default); by contrast, 
# if we gave the full path to a local pdb file, then this second argument should be online = FALSE. 
# (3) The third argument, which is named threshold, refers to the maximum distance (in Angstroms) 
# between a sulfur atom (SD) and an aromatic ring to be considered as an S-aromatic motif. For each
# methionine, the script computes the distances between its S atom and all the aromatic rings from 
# each Tyr (cY), Phe (cF) and Trp residues (The indole ring of Trp is treated as two separate rings 
# of 5 (cW1) and 6 (cW2) atoms. These distances are collected into the R matrices d_SD_Y, d_SD_F, 
# d_SD_W1 and d_SD_W2, respectively. The rownames of this matrices correspond to the methionine 
# residue and its chain, while the colnames identified the aromatic residue and its chain. 
# To summarize these results, for each methionie we identify the closest Tyr, Phe, and Trp residues 
# and the involved distances. Eventually, the program returns a dataframe where each row corresponds 
# to a methionine residue from the protein under analysis. For each methionine residue the following 
# variables (among others) are determined:
#                                    
# Met: residue number for the given methionine. 
# Met_Chain: the chain ID at which the methionine residue belongs.
# Tyr: number of the closest (to the methionine) tyrosine residue.
# Tyr_Chain: the chain ID at which that tyrosine residue belongs.
# Yd: distance (in Angstroms) between the SD and the closest tyrosine.
# Remark: diferent to "none" when more than one tyrosine are equally close.
# NumberYclosertoThreshold: as the variable's named indicates.
# The very same structure of variables is repeted for phenylalanine, the 5 atoms tryptophan ring (Trp or W1) 
# and the 6 atoms tryptophan ring (Trp.1 or W2). Beawere that for some proteins, tryptophan (e.g. 1DCY) 
# or even tyrosine may be absent.
# numberBonds: number of S-aromatic motifs formed.
# closestAromaticAt: distance (in Angstroms) between the SD and the closest aromatic ring.

SaromaticDistances <- function(pdb, online = TRUE, threshold = 7){ 
  # pdb is either de ID of a protein structure to be reached on-line or the path to a local pdb file.
  # mode is either "online" or "local" depending on where the pdb file to be used is.
  # threshold is the maximun distance in Angstroms the S atom and the aromatic ring to be considered as an S-aromatic motif
  library(Rpdb) # required packages
  library(bio3d)
  td <- threshold
  
  if (! online){
    x <- Rpdb::read.pdb(pdb) # Rpdb gets the pdb from the indicated directory.
  } else {
    y <- bio3d::read.pdb(pdb) # bio3d fetchs the structure online
    bio3d::write.pdb(y, pdb) # and writes the pdb file in the wd.
    x <- Rpdb::read.pdb(pdb) # Rpdb gets the pdb from the wd.
    # system("cmd.exe", input = paste("rm", pdb, sep=" ")) # Once used, it is deleted.
  }
  
  # Checking that there are Met, Tyr, Phe and Trp residues in our protein
  is_there_M <- TRUE %in% (x$atoms$resname == "MET")
  is_there_Y <- TRUE %in% (x$atoms$resname == "TYR")
  is_there_F <- TRUE %in% (x$atoms$resname == "PHE")
  is_there_W <- TRUE %in% (x$atoms$resname == "TRP")
  
  if (!is_there_M){ # if there is not Met the script stop with a warning
    stop("THIS PROTEIN DOES NOT CONTAIN METHIONINE")
  }
  ########## Methionine Residues
  is.Msulfur <- x$atoms$resname == "MET" & x$atoms$elename == "SD"  # it marks as TRUE the sulfur atoms from methionines
  Msulfur_subset <- x$atoms[is.Msulfur,] # Subset of the pdb object restricted to sulfur atoms from Met
  Msulfur_subset$residue <- paste(Msulfur_subset$resid, 
                                  Msulfur_subset$chainid, sep=" ") # It adds the position and chain to the data frame
  centre_Msulfur <- centres(Msulfur_subset, factor = Msulfur_subset$residue) # Object of class'coords'
  centre_Msulfur$id <- row.names(centre_Msulfur)
  centre_Msulfur$point <- "SD"
  number_M_residues <- nrow(centre_Msulfur)
  
  ########## Tyrosines Residues
  if (is_there_Y) { 
    is.Yring <- x$atoms$resname == "TYR" & x$atoms$elename %in% c("CG", "CD1",
                                                                  "CD2", "CE1", "CE2", "CZ")
    Yring_subset <- x$atoms[is.Yring,]
    Yring_subset$residue <- paste(Yring_subset$resid, Yring_subset$chainid, 
                                  sep=" ")
    # Computes the coordinates for the ring's centre:
    centre_Yring <- centres(Yring_subset, factor = Yring_subset$residue) 
    # The residue involved is identified by this variable:
    centre_Yring$id <- row.names(centre_Yring)  
    # The type of spatial point is identified (cY: centre from tyr ring):
    centre_Yring$point <- "cY" 
    number_Y_residues <- nrow(centre_Yring)
  } else { # An empty object of class "coords" is created
    centre_Yring <- centre_Msulfur[1,] # Because it need to be of class "coords"
    centre_Yring[1,] <- NA # but it should be empty: with NA 
    centre_Yring[1,5] <- "cY"
  }
  ########## Phenilalanine Residues
  if (is_there_F) { 
    is.Fring <- x$atoms$resname == "PHE" & x$atoms$elename %in% c("CG", "CD1",
                                                                  "CD2", "CE1", "CE2", "CZ")
    Fring_subset <- x$atoms[is.Fring,]
    Fring_subset$residue <- paste(Fring_subset$resid, Fring_subset$chainid, 
                                  sep=" ")
    # Computes the coordinates for the ring's centre:
    centre_Fring <- centres(Fring_subset, factor = Fring_subset$residue)
    centre_Fring$id <- row.names(centre_Fring)
    # the type of spatial point is identified (cF: centre from phe ring):
    centre_Fring$point <- "cF" 
    number_F_residues <- nrow(centre_Fring)
  }  else { # An empty object of class "coords" is created
    centre_Fring <- centre_Msulfur[1,] # Because it need to be of class "coords"
    centre_Fring[1,] <- NA # but it should be empty: with NA 
    centre_Fring[1,5] <- "cF"
  }
  ########## Tryptophan Residues
  if (is_there_W) { 
    is.Wring1 <- x$atoms$resname == "TRP" & x$atoms$elename %in% c("CG", 
                                                                   "CD1", "CD2", "CE2", "NE1")
    Wring_subset1 <- x$atoms[is.Wring1,]
    Wring_subset1$residue <- paste(Wring_subset1$resid, Wring_subset1$chainid,
                                   sep=" ")
    # Computes the coordinates for the ring's centre:
    centre_Wring1 <- centres(Wring_subset1, factor = Wring_subset1$residue)
    centre_Wring1$id <- row.names(centre_Wring1)
    number_W_residues <- nrow(centre_Wring1)
    centre_Wring1$point <- "cW1" # (cW1: centre from the 5 atoms W ring)
    
    
    is.Wring2 <- x$atoms$resname == "TRP" & x$atoms$elename %in% c("CD2", 
                                                                   "CE2", "CE3", "CZ3", "CH2", "CZ2")
    Wring_subset2 <- x$atoms[is.Wring2,]
    Wring_subset2$residue <- paste(Wring_subset2$resid, Wring_subset2$chainid,
                                   sep=" ")
    centre_Wring2 <- centres(Wring_subset2, factor = Wring_subset2$residue)
    centre_Wring2$id <- row.names(centre_Wring2)
    centre_Wring2$point <- "cW2" # (cW2: centre from the 6 atoms W ring)
  } else { # An empty object of class "coords" is created
    centre_Wring1 <- centre_Msulfur[1,] # Because this object must be of class "coord" 
    centre_Wring1[1,] <- NA # But it should be empty: NA
    centre_Wring1[1,5] <- "cW1"
    centre_Wring2 <- centre_Msulfur[1,]
    centre_Wring2[1,] <- NA 
    centre_Wring2[1,5] <- "cW2"
  }
  ########## Merging Relevant Residues
  Relevant_Residues1 <- merge(centre_Yring, centre_Fring)  
  Relevant_Residues2 <- merge(centre_Wring1, centre_Wring2)
  Relevant_Residues <- merge(Relevant_Residues1, Relevant_Residues2)
  Relevant_Residues <- merge(Relevant_Residues, centre_Msulfur)
  rm(Relevant_Residues1, Relevant_Residues2)
  
  ############# Computing Distances
  
  if (is_there_Y) { # Tyrosine
    d_vector_SD_cY <- distances(Relevant_Residues, 
                                Relevant_Residues$point == "SD", 
                                Relevant_Residues$point == "cY")
    d_SD_cY <- round(norm(d_vector_SD_cY), 1)
    rownames(d_SD_cY) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
                                           4] # To id the Met 
    colnames(d_SD_cY) <- Relevant_Residues[which(Relevant_Residues$point == "cY"),
                                           4] # To id the Tyr
  }
  
  if (is_there_F){ # Phenylalanine
    d_vector_SD_cF <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
                                Relevant_Residues$point == "cF")
    d_SD_cF <- round(norm(d_vector_SD_cF), 1)
    rownames(d_SD_cF) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
                                           4] # To id the Met 
    colnames(d_SD_cF) <- Relevant_Residues[which(Relevant_Residues$point == "cF"),
                                           4] # To id the Phe
  }
  
  if (is_there_W) { # Tryptophan
    d_vector_SD_cW1 <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
                                 Relevant_Residues$point == "cW1")
    d_SD_cW1 <- round(norm(d_vector_SD_cW1),1)
    rownames(d_SD_cW1) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
                                            4] # To id the Met 
    colnames(d_SD_cW1) <- Relevant_Residues[which(Relevant_Residues$point == "cW1"),
                                            4] # To id the Trp
    
    d_vector_SD_cW2 <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
                                 Relevant_Residues$point == "cW2")
    d_SD_cW2 <- round(norm(d_vector_SD_cW2),1)
    rownames(d_SD_cW2) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
                                            4] # To id the Met 
    colnames(d_SD_cW2) <- Relevant_Residues[which(Relevant_Residues$point == "cW2"),
                                            4] # To id the Trp
  }
  
  ########### The closest tyrosine to my methionine is at ...
  if (is_there_Y){ # When the protein contains tyrosine residues
    closestY <- data.frame(Met=numeric(), Met_Chain=factor(), Tyr=numeric(), 
                           Tyr_Chain=factor(), Yd=numeric(), Remark=character(),
                           NumberYclosertoThreshold=numeric())                           
    for (i in 1:nrow(d_SD_cY)){ # number of methionines
      met_chain <- strsplit(rownames(d_SD_cY)[i], split = " ")[[1]]
      lowerthan <- which(d_SD_cY[i,] <= td) # column indices that fulfils the cond.
      # column indice (column name gives the tyr residue and chain):
      tY <- which(d_SD_cY[i,] == min(d_SD_cY[i,])) 
      targetY <- colnames(d_SD_cY)[tY]  
      if (length(targetY) > 1){ # It may happen that a S atom is equidistant 
        # to two or more Y residues
        # It warns than more than one residue are equally close:
        v_remark = paste(targetY, collapse=";") 
        targetY <- targetY[1] # Arbitrarily, we take the first in the list
      }
      else { v_remark = "none"}
      tyr_chain <- strsplit(targetY, split = " ")[[1]] 
      closestY <- rbind(closestY, data.frame(Met=as.numeric(met_chain[1]), 
                                             Met_Chain=met_chain[2], Tyr=as.numeric(tyr_chain[1]), 
                                             Tyr_Chain=tyr_chain[2], Yd=d_SD_cY[i,targetY],
                                             Remark=v_remark, NumberYclosertoThreshold = length(lowerthan)))  
    } 
  } else { # If the protein doesn't contain tyrosine residues
    
    closestY <- data.frame(Met=numeric(), Met_Chain=factor(), Tyr=numeric(), 
                           Tyr_Chain=factor(), Yd=numeric(),Remark=character(),
                           NumberYclosertoThreshold=numeric()) 
    for (i in 1:nrow(centre_Msulfur)){ # number of Met
      met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines 
      
      closestY <- rbind(closestY, data.frame(Met=as.numeric(met_chain[1]), 
                                             Met_Chain=met_chain[2], Tyr="any", 
                                             Tyr_Chain="any", Yd=Inf, 
                                             Remark="any", NumberYclosertoThreshold = 0))
    }
  }
  ########## The closest phenylalanine to my methionine is at ...
  if (is_there_F){ 
    closestF <- data.frame(Met=numeric(), Met_Chain=factor(), Phe=numeric(), 
                           Phe_Chain=factor(), Fd=numeric(), Remark=character(),
                           NumberFclosertoThreshold=numeric())                           
    for (i in 1:nrow(d_SD_cF)){
      met_chain <- strsplit(rownames(d_SD_cF)[i], split = " ")[[1]]
      lowerthan <- which(d_SD_cF[i,] <= td) # column indices that fulfils the cond.
      # Column indice (column name gives the phe residue and chain):
      tF <- which(d_SD_cF[i,] == min(d_SD_cF[i,])) 
      targetF <- colnames(d_SD_cF)[tF]
      if (length(targetF) > 1){ # It may happen that a S atom is equidistant 
        # to two or more F residues
        # It warns that more than one residue are equally close:
        v_remark = paste(targetF, collapse=";") 
        targetF <- targetF[1] # Arbitrarily, we take the first in the list
      }
      else { v_remark = "none"}
      phe_chain <- strsplit(targetF, split = " ")[[1]] 
      closestF <- rbind(closestF, data.frame(Met=as.numeric(met_chain[1]), 
                                             Met_Chain=met_chain[2], Phe=as.numeric(phe_chain[1]), 
                                             Phe_Chain=phe_chain[2], Fd=d_SD_cF[i,targetF], 
                                             Remark=v_remark, NumberFclosertoThreshold = length(lowerthan)))  
    }
  } else { # If the protein doesn't contain phenylalanine
    
    closestF <- data.frame(Met=numeric(), Met_Chain=factor(), Phe=numeric(), 
                           Phe_Chain=factor(), Fd=numeric(),Remark=character(),
                           NumberFclosertoThreshold=numeric()) 
    for (i in 1:nrow(centre_Msulfur)){ # number of Met
      met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines   
      
      closestF <- rbind(closestF, data.frame(Met=as.numeric(met_chain[1]), 
                                             Met_Chain=met_chain[2], Phe="any", 
                                             Phe_Chain="any", Fd=Inf, 
                                             Remark="any", NumberFclosertoThreshold = 0))
    }
  }
  ########## The closest Tryptophan to my methionine is at ...
  if (is_there_W) { # If the protein contains tryptophan
    closestW1 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
                            Trp_Chain=factor(), W1d=numeric(), Remark=character(),
                            NumberW1closertoThreshold=numeric())                           
    for (i in 1:nrow(d_SD_cW1)){
      met_chain <- strsplit(rownames(d_SD_cW1)[i], split = " ")[[1]]
      lowerthan <- which(d_SD_cW1[i,] <= td) # column indices that fulfils the cond.
      # Column indice (column name gives the trp residue and chain):
      tW1 <- which(d_SD_cW1[i,] == min(d_SD_cW1[i,]))
      targetW1 <- colnames(d_SD_cW1)[tW1]      
      if (length(targetW1) > 1){ # It may happen that a S atom is equidistant
        # to two or more W residues
        # It warns that more than one residue are equally close:
        v_remark = paste(targetW1, collapse=";") 
        targetW1 <- targetW1[1] # Arbitrarily, we take the first in the list
      }
      else { v_remark = "none"}
      trp_chain <- strsplit(targetW1, split = " ")[[1]] 
      closestW1 <- rbind(closestW1, data.frame(Met=as.numeric(met_chain[1]),
                                               Met_Chain=met_chain[2], Trp=as.numeric(trp_chain[1]), 
                                               Trp_Chain=trp_chain[2], W1d=d_SD_cW1[i,targetW1], 
                                               Remark=v_remark, NumberW1closertoThreshold = length(lowerthan)))  
    }
    
    closestW2 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
                            Trp_Chain=factor(), W2d=numeric(),Remark=character(),
                            NumberW2closertoThreshold=numeric())                           
    for (i in 1:nrow(d_SD_cW2)){
      met_chain <- strsplit(rownames(d_SD_cW2)[i], split = " ")[[1]]
      lowerthan <- which(d_SD_cW2[i,] <= td) # column indices that fulfils the cond.
      # Column indice, colnames gives the trp residue (6 atoms ring) and chain: 
      tW2 <- which(d_SD_cW2[i,] == min(d_SD_cW2[i,]))
      targetW2 <- colnames(d_SD_cW2)[tW2]     
      # It may happen that a S atom is equidistant to two or more W residues:
      if (length(targetW2) > 1){ 
        # It warns than more than one residue are equally close:
        v_remark = paste(targetW2, collapse=";") 
        targetW2 <- targetW2[1] # Arbitrarily, we take the first in the list
      }
      else { v_remark = "none"}
      trp_chain <- strsplit(targetW2, split = " ")[[1]] 
      closestW2 <- rbind(closestW2, data.frame(Met=as.numeric(met_chain[1]), 
                                               Met_Chain=met_chain[2], Trp=as.numeric(trp_chain[1]),
                                               Trp_Chain=trp_chain[2], W2d=d_SD_cW2[i,targetW2], 
                                               Remark=v_remark, NumberW2closertoThreshold = length(lowerthan)))  
    }
  } else { # If the protein doesn't contain tryptophan
    
    closestW1 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
                            Trp_Chain=factor(), W1d=numeric(),Remark=character(),
                            NumberW1closertoThreshold=numeric())  
    closestW2 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
                            Trp_Chain=factor(), W2d=numeric(),Remark=character(),
                            NumberW2closertoThreshold=numeric()) 
    for (i in 1:nrow(centre_Msulfur)){ # number of Met
      met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines   
      
      closestW1 <- rbind(closestW1, data.frame(Met=as.numeric(met_chain[1]), 
                                               Met_Chain=met_chain[2], Trp="any", 
                                               Trp_Chain="any", W1d=Inf, 
                                               Remark="any", NumberW1closertoThreshold = 0))
      closestW2 <- rbind(closestW2, data.frame(Met=as.numeric(met_chain[1]), 
                                               Met_Chain=met_chain[2], Trp="any", 
                                               Trp_Chain="any", W2d=Inf, 
                                               Remark="any", NumberW2closertoThreshold = 0))
    }
  }
  ########## The closest aromatic rings to my methionine are ...
  # cA: closest-Aromatic
  cA <- data.frame(closestY, closestF, closestW1, closestW2) 
  checking <- cA$Met == cA$Met.1 & cA$Met == cA$Met.2 & cA$Met == cA$Met.3
  
  if (FALSE %in% checking) {
    stop("PROBLEMS WITH THE METHIONINE IDENTITIES WHEN MARGING RESULTS")
  } else {
    cA$Met.1 <- NULL
    cA$Met.2 <- NULL
    cA$Met.3 <- NULL
  }
  checking <- cA$Met_Chain == cA$Met_Chain.1 & cA$Met_Chain == cA$Met_Chain.2 &
    cA$Met_Chain == cA$Met_Chain.3
  
  if (FALSE %in% checking) {
    stop("PROBLEMS WITH THE CHAIN IDENTITIES WHEN MARGING RESULTS")
  } else {
    cA$Met_Chain.1 <- NULL
    cA$Met_Chain.2 <- NULL
    cA$Met_Chain.3 <- NULL
  }
  
  cA$numberBonds <- -999 # by default. It will be changed by the program in due course
  cA$closestAromaticAt <- -999 
  cA$Wd <- -999
  for (i in 1:nrow(cA)){
    NumberW <- max(cA$NumberW1closertoThreshold[i], cA$NumberW2closertoThreshold[i])
    cA$numberBonds[i] <- cA$NumberYclosertoThreshold[i]+ cA$NumberFclosertoThreshold[i] + NumberW
    cA$closestAromaticAt[i] <- min(cA$Yd[i],cA$Fd[i],cA$W1d[i],cA$W2d[i])
    cA$Wd[i] <- min(cA$W1d[i], cA$W2d[i])
  }  
  
  cA <- cA[order(cA[,1]),] # Order the rows by increasing Met position
  return(cA)
}
#################################################################################################
# MF <- function(pdb, online = TRUE, threshold = 7, res_aro = F){ 
#   # pdb is either de ID of a protein structure to be reached on-line or the path to a local pdb file.
#   # mode is either "online" or "local" depending on where the pdb file to be used is.
#   # threshold is the maximun distance in Angstroms the S atom and the aromatic ring to be considered as an S-aromatic motif
#   library(Rpdb) # required packages
#   library(bio3d)
#   td <- threshold
#   
#   if (! online){
#     x <- Rpdb::read.pdb(pdb) # Rpdb gets the pdb from the indicated directory.
#   } else {
#     y <- bio3d::read.pdb(pdb) # bio3d fetchs the structure online
#     bio3d::write.pdb(y, pdb) # and writes the pdb file in the wd.
#     x <- Rpdb::read.pdb(pdb) # Rpdb gets the pdb from the wd.
#     # system("cmd.exe", input = paste("rm", pdb, sep=" ")) # Once used, it is deleted.
#   }
#   
#   # Checking that there are Met, Tyr, Phe and Trp residues in our protein
#   is_there_M <- TRUE %in% (x$atoms$resname == "MET")
#   is_there_Y <- TRUE %in% (x$atoms$resname == "TYR")
#   is_there_F <- TRUE %in% (x$atoms$resname == "PHE")
#   is_there_W <- TRUE %in% (x$atoms$resname == "TRP")
#   
#   if (!is_there_M){ # if there is not Met the script stop with a warning
#     stop("THIS PROTEIN DOES NOT CONTAIN METHIONINE")
#   }
#   ########## Methionine Residues
#   is.Msulfur <- x$atoms$resname == "MET" & x$atoms$elename == "SD"  # it marks as TRUE the sulfur atoms from methionines
#   Msulfur_subset <- x$atoms[is.Msulfur,] # Subset of the pdb object restricted to sulfur atoms from Met
#   Msulfur_subset$residue <- paste(Msulfur_subset$resid, 
#                                   Msulfur_subset$chainid, sep=" ") # It adds the position and chain to the data frame
#   centre_Msulfur <- centres(Msulfur_subset, factor = Msulfur_subset$residue) # Object of class'coords'
#   centre_Msulfur$id <- row.names(centre_Msulfur)
#   centre_Msulfur$point <- "SD"
#   number_M_residues <- nrow(centre_Msulfur)
#   
#   ########## Tyrosines Residues
#   if (is_there_Y) { 
#     is.Yring <- x$atoms$resname == "TYR" & x$atoms$elename %in% c("CG", "CD1",
#                                                                   "CD2", "CE1", "CE2", "CZ")
#     Yring_subset <- x$atoms[is.Yring,]
#     Yring_subset$residue <- paste(Yring_subset$resid, Yring_subset$chainid, 
#                                   sep=" ")
#     # Computes the coordinates for the ring's centre:
#     centre_Yring <- centres(Yring_subset, factor = Yring_subset$residue) 
#     # The residue involved is identified by this variable:
#     centre_Yring$id <- row.names(centre_Yring)  
#     # The type of spatial point is identified (cY: centre from tyr ring):
#     centre_Yring$point <- "cY" 
#     number_Y_residues <- nrow(centre_Yring)
#   } else { # An empty object of class "coords" is created
#     centre_Yring <- centre_Msulfur[1,] # Because it need to be of class "coords"
#     centre_Yring[1,] <- NA # but it should be empty: with NA 
#     centre_Yring[1,5] <- "cY"
#   }
#   ########## Phenilalanine Residues
#   if (is_there_F) { 
#     is.Fring <- x$atoms$resname == "PHE" & x$atoms$elename %in% c("CG", "CD1",
#                                                                   "CD2", "CE1", "CE2", "CZ")
#     Fring_subset <- x$atoms[is.Fring,]
#     Fring_subset$residue <- paste(Fring_subset$resid, Fring_subset$chainid, 
#                                   sep=" ")
#     # Computes the coordinates for the ring's centre:
#     centre_Fring <- centres(Fring_subset, factor = Fring_subset$residue)
#     centre_Fring$id <- row.names(centre_Fring)
#     # the type of spatial point is identified (cF: centre from phe ring):
#     centre_Fring$point <- "cF" 
#     number_F_residues <- nrow(centre_Fring)
#   }  else { # An empty object of class "coords" is created
#     centre_Fring <- centre_Msulfur[1,] # Because it need to be of class "coords"
#     centre_Fring[1,] <- NA # but it should be empty: with NA 
#     centre_Fring[1,5] <- "cF"
#   }
#   ########## Tryptophan Residues
#   if (is_there_W) { 
#     is.Wring1 <- x$atoms$resname == "TRP" & x$atoms$elename %in% c("CG", 
#                                                                    "CD1", "CD2", "CE2", "NE1")
#     Wring_subset1 <- x$atoms[is.Wring1,]
#     Wring_subset1$residue <- paste(Wring_subset1$resid, Wring_subset1$chainid,
#                                    sep=" ")
#     # Computes the coordinates for the ring's centre:
#     centre_Wring1 <- centres(Wring_subset1, factor = Wring_subset1$residue)
#     centre_Wring1$id <- row.names(centre_Wring1)
#     number_W_residues <- nrow(centre_Wring1)
#     centre_Wring1$point <- "cW1" # (cW1: centre from the 5 atoms W ring)
#     
#     
#     is.Wring2 <- x$atoms$resname == "TRP" & x$atoms$elename %in% c("CD2", 
#                                                                    "CE2", "CE3", "CZ3", "CH2", "CZ2")
#     Wring_subset2 <- x$atoms[is.Wring2,]
#     Wring_subset2$residue <- paste(Wring_subset2$resid, Wring_subset2$chainid,
#                                    sep=" ")
#     centre_Wring2 <- centres(Wring_subset2, factor = Wring_subset2$residue)
#     centre_Wring2$id <- row.names(centre_Wring2)
#     centre_Wring2$point <- "cW2" # (cW2: centre from the 6 atoms W ring)
#   } else { # An empty object of class "coords" is created
#     centre_Wring1 <- centre_Msulfur[1,] # Because this object must be of class "coord" 
#     centre_Wring1[1,] <- NA # But it should be empty: NA
#     centre_Wring1[1,5] <- "cW1"
#     centre_Wring2 <- centre_Msulfur[1,]
#     centre_Wring2[1,] <- NA 
#     centre_Wring2[1,5] <- "cW2"
#   }
#   ########## Merging Relevant Residues
#   Relevant_Residues1 <- merge(centre_Yring, centre_Fring)  
#   Relevant_Residues2 <- merge(centre_Wring1, centre_Wring2)
#   Relevant_Residues <- merge(Relevant_Residues1, Relevant_Residues2)
#   Relevant_Residues <- merge(Relevant_Residues, centre_Msulfur)
#   rm(Relevant_Residues1, Relevant_Residues2)
#   
#   ############# Computing Distances
#   
#   if (is_there_Y) { # Tyrosine
#     d_vector_SD_cY <- distances(Relevant_Residues, 
#                                 Relevant_Residues$point == "SD", 
#                                 Relevant_Residues$point == "cY")
#     d_SD_cY <- round(norm(d_vector_SD_cY), 1)
#     rownames(d_SD_cY) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
#                                            4] # To id the Met 
#     colnames(d_SD_cY) <- Relevant_Residues[which(Relevant_Residues$point == "cY"),
#                                            4] # To id the Tyr
#   }
#   
#   if (is_there_F){ # Phenylalanine
#     d_vector_SD_cF <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
#                                 Relevant_Residues$point == "cF")
#     d_SD_cF <- round(norm(d_vector_SD_cF), 1)
#     rownames(d_SD_cF) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
#                                            4] # To id the Met 
#     colnames(d_SD_cF) <- Relevant_Residues[which(Relevant_Residues$point == "cF"),
#                                            4] # To id the Phe
#   }
#   
#   if (is_there_W) { # Tryptophan
#     d_vector_SD_cW1 <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
#                                  Relevant_Residues$point == "cW1")
#     d_SD_cW1 <- round(norm(d_vector_SD_cW1),1)
#     rownames(d_SD_cW1) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
#                                             4] # To id the Met 
#     colnames(d_SD_cW1) <- Relevant_Residues[which(Relevant_Residues$point == "cW1"),
#                                             4] # To id the Trp
#     
#     d_vector_SD_cW2 <- distances(Relevant_Residues, Relevant_Residues$point == "SD", 
#                                  Relevant_Residues$point == "cW2")
#     d_SD_cW2 <- round(norm(d_vector_SD_cW2),1)
#     rownames(d_SD_cW2) <- Relevant_Residues[which(Relevant_Residues$point == "SD"),
#                                             4] # To id the Met 
#     colnames(d_SD_cW2) <- Relevant_Residues[which(Relevant_Residues$point == "cW2"),
#                                             4] # To id the Trp
#   }
#   
#   ########### The closest tyrosine to my methionine is at ...
#   if (is_there_Y){ # When the protein contains tyrosine residues
#     closestY <- data.frame(Met=numeric(), Met_Chain=factor(), Tyr=numeric(), 
#                            Tyr_Chain=factor(), Yd=numeric(), Remark=character(),
#                            NumberYclosertoThreshold=numeric())                           
#     for (i in 1:nrow(d_SD_cY)){ # number of methionines
#       met_chain <- strsplit(rownames(d_SD_cY)[i], split = " ")[[1]]
#       lowerthan <- which(d_SD_cY[i,] <= td) # column indices that fulfils the cond.
#       # column indice (column name gives the tyr residue and chain):
#       tY <- which(d_SD_cY[i,] == min(d_SD_cY[i,])) 
#       targetY <- colnames(d_SD_cY)[tY]  
#       if (length(targetY) > 1){ # It may happen that a S atom is equidistant 
#         # to two or more Y residues
#         # It warns than more than one residue are equally close:
#         v_remark = paste(targetY, collapse=";") 
#         targetY <- targetY[1] # Arbitrarily, we take the first in the list
#       }
#       else { v_remark = "none"}
#       tyr_chain <- strsplit(targetY, split = " ")[[1]] 
#       closestY <- rbind(closestY, data.frame(Met=as.numeric(met_chain[1]), 
#                                              Met_Chain=met_chain[2], Tyr=as.numeric(tyr_chain[1]), 
#                                              Tyr_Chain=tyr_chain[2], Yd=d_SD_cY[i,targetY],
#                                              Remark=v_remark, NumberYclosertoThreshold = length(lowerthan)))  
#     } 
#   } else { # If the protein doesn't contain tyrosine residues
#     
#     closestY <- data.frame(Met=numeric(), Met_Chain=factor(), Tyr=numeric(), 
#                            Tyr_Chain=factor(), Yd=numeric(),Remark=character(),
#                            NumberYclosertoThreshold=numeric()) 
#     for (i in 1:nrow(centre_Msulfur)){ # number of Met
#       met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines 
#       
#       closestY <- rbind(closestY, data.frame(Met=as.numeric(met_chain[1]), 
#                                              Met_Chain=met_chain[2], Tyr="any", 
#                                              Tyr_Chain="any", Yd=Inf, 
#                                              Remark="any", NumberYclosertoThreshold = 0))
#     }
#   }
#   ########## The closest phenylalanine to my methionine is at ...
#   if (is_there_F){ 
#     closestF <- data.frame(Met=numeric(), Met_Chain=factor(), Phe=numeric(), 
#                            Phe_Chain=factor(), Fd=numeric(), Remark=character(),
#                            NumberFclosertoThreshold=numeric())                           
#     for (i in 1:nrow(d_SD_cF)){
#       met_chain <- strsplit(rownames(d_SD_cF)[i], split = " ")[[1]]
#       lowerthan <- which(d_SD_cF[i,] <= td) # column indices that fulfils the cond.
#       # Column indice (column name gives the phe residue and chain):
#       tF <- which(d_SD_cF[i,] == min(d_SD_cF[i,])) 
#       targetF <- colnames(d_SD_cF)[tF]
#       if (length(targetF) > 1){ # It may happen that a S atom is equidistant 
#         # to two or more F residues
#         # It warns that more than one residue are equally close:
#         v_remark = paste(targetF, collapse=";") 
#         targetF <- targetF[1] # Arbitrarily, we take the first in the list
#       }
#       else { v_remark = "none"}
#       phe_chain <- strsplit(targetF, split = " ")[[1]] 
#       closestF <- rbind(closestF, data.frame(Met=as.numeric(met_chain[1]), 
#                                              Met_Chain=met_chain[2], Phe=as.numeric(phe_chain[1]), 
#                                              Phe_Chain=phe_chain[2], Fd=d_SD_cF[i,targetF], 
#                                              Remark=v_remark, NumberFclosertoThreshold = length(lowerthan)))  
#     }
#   } else { # If the protein doesn't contain phenylalanine
#     
#     closestF <- data.frame(Met=numeric(), Met_Chain=factor(), Phe=numeric(), 
#                            Phe_Chain=factor(), Fd=numeric(),Remark=character(),
#                            NumberFclosertoThreshold=numeric()) 
#     for (i in 1:nrow(centre_Msulfur)){ # number of Met
#       met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines   
#       
#       closestF <- rbind(closestF, data.frame(Met=as.numeric(met_chain[1]), 
#                                              Met_Chain=met_chain[2], Phe="any", 
#                                              Phe_Chain="any", Fd=Inf, 
#                                              Remark="any", NumberFclosertoThreshold = 0))
#     }
#   }
#   ########## The closest Tryptophan to my methionine is at ...
#   if (is_there_W) { # If the protein contains tryptophan
#     closestW1 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
#                             Trp_Chain=factor(), W1d=numeric(), Remark=character(),
#                             NumberW1closertoThreshold=numeric())                           
#     for (i in 1:nrow(d_SD_cW1)){
#       met_chain <- strsplit(rownames(d_SD_cW1)[i], split = " ")[[1]]
#       lowerthan <- which(d_SD_cW1[i,] <= td) # column indices that fulfils the cond.
#       # Column indice (column name gives the trp residue and chain):
#       tW1 <- which(d_SD_cW1[i,] == min(d_SD_cW1[i,]))
#       targetW1 <- colnames(d_SD_cW1)[tW1]      
#       if (length(targetW1) > 1){ # It may happen that a S atom is equidistant
#         # to two or more W residues
#         # It warns that more than one residue are equally close:
#         v_remark = paste(targetW1, collapse=";") 
#         targetW1 <- targetW1[1] # Arbitrarily, we take the first in the list
#       }
#       else { v_remark = "none"}
#       trp_chain <- strsplit(targetW1, split = " ")[[1]] 
#       closestW1 <- rbind(closestW1, data.frame(Met=as.numeric(met_chain[1]),
#                                                Met_Chain=met_chain[2], Trp=as.numeric(trp_chain[1]), 
#                                                Trp_Chain=trp_chain[2], W1d=d_SD_cW1[i,targetW1], 
#                                                Remark=v_remark, NumberW1closertoThreshold = length(lowerthan)))  
#     }
#     
#     closestW2 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
#                             Trp_Chain=factor(), W2d=numeric(),Remark=character(),
#                             NumberW2closertoThreshold=numeric())                           
#     for (i in 1:nrow(d_SD_cW2)){
#       met_chain <- strsplit(rownames(d_SD_cW2)[i], split = " ")[[1]]
#       lowerthan <- which(d_SD_cW2[i,] <= td) # column indices that fulfils the cond.
#       # Column indice, colnames gives the trp residue (6 atoms ring) and chain: 
#       tW2 <- which(d_SD_cW2[i,] == min(d_SD_cW2[i,]))
#       targetW2 <- colnames(d_SD_cW2)[tW2]     
#       # It may happen that a S atom is equidistant to two or more W residues:
#       if (length(targetW2) > 1){ 
#         # It warns than more than one residue are equally close:
#         v_remark = paste(targetW2, collapse=";") 
#         targetW2 <- targetW2[1] # Arbitrarily, we take the first in the list
#       }
#       else { v_remark = "none"}
#       trp_chain <- strsplit(targetW2, split = " ")[[1]] 
#       closestW2 <- rbind(closestW2, data.frame(Met=as.numeric(met_chain[1]), 
#                                                Met_Chain=met_chain[2], Trp=as.numeric(trp_chain[1]),
#                                                Trp_Chain=trp_chain[2], W2d=d_SD_cW2[i,targetW2], 
#                                                Remark=v_remark, NumberW2closertoThreshold = length(lowerthan)))  
#     }
#   } else { # If the protein doesn't contain tryptophan
#     
#     closestW1 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
#                             Trp_Chain=factor(), W1d=numeric(),Remark=character(),
#                             NumberW1closertoThreshold=numeric())  
#     closestW2 <- data.frame(Met=numeric(), Met_Chain=factor(), Trp=numeric(), 
#                             Trp_Chain=factor(), W2d=numeric(),Remark=character(),
#                             NumberW2closertoThreshold=numeric()) 
#     for (i in 1:nrow(centre_Msulfur)){ # number of Met
#       met_chain <- strsplit(centre_Msulfur$id[i], split = " ")[[1]] # position-chain of the methionines   
#       
#       closestW1 <- rbind(closestW1, data.frame(Met=as.numeric(met_chain[1]), 
#                                                Met_Chain=met_chain[2], Trp="any", 
#                                                Trp_Chain="any", W1d=Inf, 
#                                                Remark="any", NumberW1closertoThreshold = 0))
#       closestW2 <- rbind(closestW2, data.frame(Met=as.numeric(met_chain[1]), 
#                                                Met_Chain=met_chain[2], Trp="any", 
#                                                Trp_Chain="any", W2d=Inf, 
#                                                Remark="any", NumberW2closertoThreshold = 0))
#     }
#   }
#   ########## The closest aromatic rings to my methionine are ...
#   # cA: closest-Aromatic
#   cA <- data.frame(closestY, closestF, closestW1, closestW2) 
#   checking <- cA$Met == cA$Met.1 & cA$Met == cA$Met.2 & cA$Met == cA$Met.3
#   
#   if (FALSE %in% checking) {
#     stop("PROBLEMS WITH THE METHIONINE IDENTITIES WHEN MARGING RESULTS")
#   } else {
#     cA$Met.1 <- NULL
#     cA$Met.2 <- NULL
#     cA$Met.3 <- NULL
#   }
#   checking <- cA$Met_Chain == cA$Met_Chain.1 & cA$Met_Chain == cA$Met_Chain.2 &
#     cA$Met_Chain == cA$Met_Chain.3
#   
#   if (FALSE %in% checking) {
#     stop("PROBLEMS WITH THE CHAIN IDENTITIES WHEN MARGING RESULTS")
#   } else {
#     cA$Met_Chain.1 <- NULL
#     cA$Met_Chain.2 <- NULL
#     cA$Met_Chain.3 <- NULL
#   }
#   
#   cA$numberBonds <- -999 # by default. It will be changed by the program in due course
#   cA$closestAromaticAt <- -999 
#   cA$Wd <- -999
#   for (i in 1:nrow(cA)){
#     NumberW <- max(cA$NumberW1closertoThreshold[i], cA$NumberW2closertoThreshold[i])
#     cA$numberBonds[i] <- cA$NumberYclosertoThreshold[i]+ cA$NumberFclosertoThreshold[i] + NumberW
#     cA$closestAromaticAt[i] <- min(cA$Yd[i],cA$Fd[i],cA$W1d[i],cA$W2d[i])
#     cA$Wd[i] <- min(cA$W1d[i], cA$W2d[i])
#   }  
#   
#   cA <- cA[order(cA[,1]),] # Order the rows by increasing Met position
#   if (res_aro == T){
#     if (!is_there_F){
#       stop("PROTEIN WITHOUT PHENILALANINE RESIDUES")
#     } else {
#       d_SD <- d_SD_cF
#     }
#     if (is_there_Y){
#       d_SD <- cbind(d_SD, d_SD_cY)
#     }
#     if (is_there_W){
#       d_SD <- cbind(d_SD_cW1, d_SD_cW2)
#     }
#     return(d_SD)
#   } else {
#     return(cA)
#   }
#   
# }
# #################################################################################################

#############     PHOSPHOSITES ANALYSIS     ######################################################
# This script takes as argument the PDB identifier (PDB) 
# of the protein under study. The program identifies the methionines present in the 
# current structure, and computes the distance in Angstroms between the sulfur atom from 
# the methionine and the oxygen atom from the phosphorylatable residue, for all the 
# phosphosites detected in the protein. This function returns a dataframe with as many
# rows as methionine are found in the protein. Each column gives the following information: 

# ptmDF$closest.ptm: distance in A at the closest PTM site.
# ptmDF$closer10A: number of  PTM sites within a sphere of radio 10 A.
# ptmDF$away.ptm: distance, in number of residues, to the closest PTM site.
# ptmDF$closer10res: number of  PTM sites within an interval 21 residues centred at the methionine.
# ptmDF$closest_pS: distance, in number of residues, to the closest pSer site
# ptmDF$closest_pT: distance, in number of residues, to the closest pThr site
# ptmDF$closest_pY: distance, in number of residues, to the closest pTyr site


phosphosite <- function(pdbID, Chain, email){
  
  library(bio3d)
  # source("/Users/JCA/Dropbox/Programacion/R/R_ToolBox/pdb2acc.R")
  source("pdb2acc.R")
  ACC <- pdb2acc(pdbID, email)[1]
  # source("/Users/JCA/Dropbox/Programacion/R/R_ToolBox/dif.R")
  source("dif.R")
  dif <- dif(pdbID, Chain, email)
  # dif <- # Proteína - # PDB (when the numeration from the UniProt and PDB files are not coincident)
  
  # From UniProt:
  download.file(paste("http://www.uniprot.org/uniprot/", ACC, ".fasta", sep=""), destfile="./temp.fas")
  uniprot <- read.fasta("./temp.fas")
  system("rm temp.fas") # Unix
  # system("cmd.exe", input = "del temp.fas") # Windows
  
  # From PDB:
  mypdb <- bio3d::read.pdb(pdbID)
  pdb.seq <- aa321(mypdb$atom$resid[which(mypdb$atom$elety == "CA" & mypdb$atom$type == "ATOM")])
  pos <- mypdb$atom[which(mypdb$atom$elety == "CA"  & mypdb$atom$type == "ATOM"),]$resno # Residue numbers.
  names(pdb.seq) <- pos + dif # Si es homooligomérica la dif es aplicable y si es heterooligomérica la dif es irrelevante.
  
  pdb.chain.seq <- aa321(mypdb$atom$resid[which(mypdb$atom$elety == "CA" & mypdb$atom$type == "ATOM" & mypdb$atom$chain == Chain)])
  pos.chain <- mypdb$atom[which(mypdb$atom$elety == "CA"  & mypdb$atom$type == "ATOM" & mypdb$atom$chain == Chain),]$resno # Residue numbers. 
  names(pdb.chain.seq) <- pos.chain + dif
  
  # Methionines found in the PDB:
  Met <- names(pdb.chain.seq[which(pdb.chain.seq == "M")]) # Posiciones de las Met en la cadena relevante.
  Met <- as.numeric(Met)
  M.pdb.seq <- paste(pdb.chain.seq[which(pdb.chain.seq == "M")], Met, sep="")
  
  # From PhosphoSitePlus:
  # all <- read.table("/Users/JCA/Dropbox/Programacion/R/R_ToolBox/Phosphorylation_site_dataset", 
  #                   header =TRUE, sep="\t", skip=3)
  all <- read.table("Phosphorylation_site_dataset", 
                    header =TRUE, sep="\t", skip=3)
  hall <- all[which(all$ORGANISM == "human"),] # Keep only human proteins
  rm(all)
  ptm <- as.character(hall[which(hall$ACC_ID == ACC),]$MOD_RSD)  # Phosphosites (Residue + Position + Modification).
  
  if (length(ptm) > 0){ 
    ptm <- sapply(ptm, function(x) gsub("-p", "", x)) # Removing the end '-p' which alludes to the modification.
    ptm <- unname(ptm)
    
    ptm.res <- sapply(ptm, function(x) substr(x, start=1, stop=1)) # Residues.
    ptm.pos <- sapply(ptm, function(x) gsub(substr(x, start=1, stop=1), "", x)) # Positions.
    
    ptm.ser <- ptm[which(sapply(ptm, function(x) substr(x, start=1, stop=1)) == "S")]
    pSer.pos <- as.numeric(sapply(ptm.ser, function(x) gsub(substr(x, start=1, stop=1), "", x)))
    
    ptm.thr <- ptm[which(sapply(ptm, function(x) substr(x, start=1, stop=1)) == "T")]
    pThr.pos <- as.numeric(sapply(ptm.thr, function(x) gsub(substr(x, start=1, stop=1), "", x)))
    
    ptm.tyr <- ptm[which(sapply(ptm, function(x) substr(x, start=1, stop=1)) == "Y")]
    pTyr.pos <- as.numeric(sapply(ptm.tyr, function(x) gsub(substr(x, start=1, stop=1), "", x)))
    
    # Checking that the phosphosites appear in the structure (PDB):
    ptm_in_pdb <- pdb.seq[which(names(pdb.seq) %in% ptm.pos)]
    ptm_in_pdb <- paste(ptm_in_pdb, names(ptm_in_pdb), sep="")  
    ptm_in_pdb <- ptm_in_pdb[which(ptm_in_pdb %in% ptm)] # Make sure they are ptm when dealing with heterooligomers.
    
    ## ------- Ancillary function to compute distances between residues ----------------- ## 
    PTM_distances <- function(x, y){ # Ej. x = T1058, y = M2288
      x.res <- gsub("[0-9]+", "", x) # Residue x (Ej. T)
      y.res <- gsub("[0-9]+", "", y) # Residue y (Ej. M)
      x.pos <- as.numeric(gsub("[A-Z]+", "", x)) - dif # Position of residue x (Ej. 1058)
      y.pos <- as.numeric(gsub("[A-Z]+", "", y)) - dif # Position of residue y (Ej. 2288)
      
      modifiable.res <- c("M","S","T","Y") 
      modifiable.atom <- c("SD", "OG", "OG1", "OH")
      x.atom <- modifiable.atom[which(modifiable.res == x.res)] # Ej. OG1
      y.atom <- modifiable.atom[which(modifiable.res == y.res)] # Ej. SD 
      
      # coordinate x.atom
      ax <- mypdb$atom$x[which(mypdb$atom$resno == x.pos & mypdb$atom$elety == x.atom)]
      bx <- mypdb$atom$y[which(mypdb$atom$resno == x.pos & mypdb$atom$elety == x.atom)]
      cx <- mypdb$atom$z[which(mypdb$atom$resno == x.pos & mypdb$atom$elety == x.atom)]
      # coordinate y.atom
      ay <- mypdb$atom$x[which(mypdb$atom$resno == y.pos & mypdb$atom$elety == y.atom)]
      by <- mypdb$atom$y[which(mypdb$atom$resno == y.pos & mypdb$atom$elety == y.atom)]
      cy <- mypdb$atom$z[which(mypdb$atom$resno == y.pos & mypdb$atom$elety == y.atom)]
      # distance
      distance <- sqrt( (ax - ay)^2 + (bx - by)^2 + (cx - cy)^2 )
      
      return(distance) 
    }
    ## ----------------------------------------------------------------------------- ##
    
    ### Related to the tertiary structure: 
    ptmDF <- data.frame(Met) # Make a DF to contain the results
    ptmDF$closest.ptm <- NA # By default. It will be modified when required. # PTM closest to Met.
    ptmDF$closer10A <- NA # By default. It will be modified when required. # Number PTM sites in a sphere of radious 10 A.
    for (i in 1:length(M.pdb.seq)){
      dist_M_ptm <- c()
      for (j in 1:length(ptm_in_pdb)){
        distancia <- round(PTM_distances(M.pdb.seq[i], ptm_in_pdb[j]), 1)
        dist_M_ptm <- c(dist_M_ptm, distancia)
      }
      ptmDF$closest.ptm[i] <- min(dist_M_ptm) # PTM closest to Met.
      ptmDF$closer10A[i] <- sum(dist_M_ptm <= 10) # Number of PTM sites in a sphere of radious 10 A.
    }
    
    ### Related to the primary structure:
    ptmDF$away.ptm <- NA # By default. It will be modified when required. # Distance in residues between Met and the closest PTM site.
    ptmDF$closer10res <- NA # By default. It will be modified when required. # Number of PTM sites within the interval o radious 10 residues.
    for (i in 1:length(Met)){
      away_M_ptm <- c() # Number residues between Met and PTM.
      for (j in 1:length(ptm.pos)){
        d <- abs(as.numeric(Met[i])-as.numeric(ptm.pos[j]))
        away_M_ptm <- c(away_M_ptm, d)
      }
      ptmDF$away.ptm[i] <- min(away_M_ptm) # PTM closest, in the primary structure, to Met.
      ptmDF$closer10res[i] <- sum(away_M_ptm < 11) # Number PTM sites in the 10 residues radious interval.
    }
    
    # Closest pSer
    MetpSer <- sapply(pSer.pos, function(x) abs(x-Met))
    if (is.null(dim(MetpSer))){
      if (is.list(MetpSer)) MetpSer <- unlist(MetpSer)
      ptmDF$closest_pS <- min(MetpSer)
    } else {
      ptmDF$closest_pS <- apply(MetpSer, 1, min)
    }
    # Closest pThr
    MetpThr <- sapply(pThr.pos, function(x) abs(x-Met))
    if (is.null(dim(MetpThr))){
      if (is.list(MetpThr)) MetpThr <- unlist(MetpThr)
      ptmDF$closest_pT <- min(MetpThr)
    } else {
      ptmDF$closest_pT <- apply(MetpThr, 1, min)
    }
    # Closest pTyr
    MetpTyr <- sapply(pTyr.pos, function(x) abs(x-Met))
    if (is.null(dim(MetpTyr))){
      if (is.list(MetpTyr)) MetpTyr <- unlist(MetpTyr)
      ptmDF$closest_pY <- min(MetpTyr)
    } else {
      ptmDF$closest_pY <- apply(MetpTyr, 1, min)
    }
    
    return(ptmDF)
  } else {
    return("non-phosphoprotein")
  } 
}
#######################################################################################################