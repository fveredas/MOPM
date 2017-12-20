predict.MetO <- function(pdbID, chain, email, cut){
  
  library(bio3d)
  source("./ToolBox.R") 
  
  mypdb <- bio3d::read.pdb(pdbID)
  PDBchain <- paste(pdbID, "_", chain, sep="") # Expression of the type ‘pdbId_chainId’.
  
  # When the protein exhibits cuaternary structure, the following command saves, in
  # the working directory, the pdb file for the whole protein (pdbID.pdb file), as well 
  # as a single pdb file just for the polypeptide chain marked as target (pdbID_chain.pdb file).
  bio3d::pdbsplit(get.pdb(pdbID), ids = PDBchain, path = "./") 
  
  # --------- Checking that the protein being analyzed contains methionine residues
  m <- mypdb$atom$resno[which(mypdb$atom$elety == "CA" & mypdb$atom$chain == chain & mypdb$atom$resid == "MET")]
  if (length(m) == 0){
    stop("THERE IS NOT METHIONINE IN THE CHOSEN PDB STRUCTURE!")
  }
  
  # -------------------- Structure of the output dataframe
  Features <- as.data.frame(matrix(rep(NA, length(m)*59), nrow = length(m)))
  names(Features) <- c("Met",
                       paste("NT_", LETTERS[c(1, 3:9, 11:14, 16:20, 22:23, 25)], sep=''),
                       paste("CT_", LETTERS[c(1, 3:9, 11:14, 16:20, 22, 25)], sep=''),
                       "entropy", "sd.entropy",
                       "SASA.pdb", "SASA.chain",
                       "nY.chain", "closestAro.chain",
                       "nF.pdb", "Fd.pdb", "nW.pdb",
                       "numberBonds.chain",
                       "Met2Y", "Met2S", "away.ptm", "closer10A.pdb",
                       "Met2S_PTM",
                       "Bfactor",
                       "dpx",
                       "prediction", "probability")
  
  # -------------------- Features related to the primary structure
  ACC <- pdb2acc(pdbID, email)[1] # Getting the UniProt accession 
  d <- dif(pdbID, chain, email)
  myseq <- get.seq(ACC)[[2]]
  system("rm seqs.fasta") # Unix
  # system("cmd.exe", input = del seqs.fasta") # Windows
  
  m <- m + d # Enumeration according to UniProt seq   
  Features$Met <- m
  
  for (i in 1:length(m)){
    NT_X <- CT_X <- c()
    
    for (aa in LETTERS[c(1, 3:9, 11:14, 16:20, 22:23, 25)]){
      x <- which(myseq == aa)
      nt <- m[i]-x
      nt <- nt[which(nt>0)]
      if (length(nt)==0){
        NT_X <- c(NT_X, NA)
      }
      else {
        NT_X <- c(NT_X, min(nt))
      }
      
      ct <- x-m[i]
      ct <- ct[which(ct>0)]
      if (length(ct)==0){
        CT_X <- c(CT_X, NA)
      }
      else {
        CT_X <- c(CT_X, min(ct))
      }
    }
    
    Features[i, 2:21] <- NT_X[1:20]
    Features[i, 22:40] <- CT_X[c(1:18, 20)] 
  }
  
  # ------------- Evolutionary features
  target <- paste(myseq, collapse="")
  
  shan8 <- eightShannon(target, ACC)
  pos_acc <- shan8[[1]] # position of methionines in the uniprot sequence
  
  if(!is.na(shan8)){
    for (i in 1:length(pos_acc)){
      for (j in 1:nrow(Features)){
        if (pos_acc[i] == Features$Met[j]){
          Features$entropy[j] <- round(shan8[[2]][i], 3)
        }
      }
    }
    Features$sd.entropy <- round(shan8[[4]], 4)
  }
  
  # ------------- Solvent Accessible Surface Area (SASA)
  a <- accDSSP(pdbID, 0, met=T, keepPDB = T) # MODIFICADO para adaptar a la nueva accDSSP función
  a <- a[which(a$Chain == chain),]
  Features$SASA.chain <- round(a$acc.subunit, 2)
  Features$SASA.pdb <- round(a$acc.complex, 2)
  rm(a)
  
  
  # -------Features realated to the S-aromatic motifs 
  d.chain <- SaromaticDistances(paste("./", PDBchain, ".pdb", sep=""), online=FALSE, threshold=7)
  d.chain$Met <- d.chain$Met + d
  for (i in 1:nrow(Features)){ # Because some sulfur atoms may not be resolved in the structure
    for (j in 1:nrow(d.chain)){
      if (Features$Met[i] == d.chain$Met[j]){
        Features$nY.chain[i] <- d.chain$NumberYclosertoThreshold[j]
        Features$closestAro.chain[i] <- d.chain$closestAromaticAt[j]
        Features$numberBonds.chain[i] <- d.chain$numberBonds[j]
      }
    }
  }
  
  d.pdb <- SaromaticDistances(paste("./", pdbID, ".pdb", sep=""), online=FALSE, threshold=7)
  d.pdb <- d.pdb[which(d.pdb$Met_Chain == chain),]
  d.pdb$Met <- d.pdb$Met + d
  for (i in 1:nrow(Features)){ # Because some sulfur atoms may not be resolved in the structure
    for (j in 1:nrow(d.pdb)){
      if (Features$Met[i] == d.pdb$Met[j]){
        Features$nF.pdb[i] <- d.pdb$NumberFclosertoThreshold[j]
        Features$Fd.pdb[i] <- d.pdb$Fd[j]
        Features$nW.pdb[i] <- d.pdb$NumberW1closertoThreshold[j]
      }
    }
  }
  rm(d.chain, d.pdb)
  
  
  # ----- Features related with the crosstalk between methionine sulfoxidation and O-phosphorylation
  for (i in 1:length(Features$Met)){
    Features$Met2Y[i] <- min(Features$NT_Y[i], Features$CT_Y[i], na.rm = TRUE)
    Features$Met2S[i] <- min(Features$NT_S[i], Features$CT_S[i], na.rm = TRUE)
  }
  p <- phosphosite(pdbID, chain, email)
  if (p != "non-phosphoprotein"){
    Features$away.ptm <- p$away.ptm
    Features$closer10A.pdb <- p$closer10A
    Features$Met2S_PTM <- p$closest_pS
    p <- "phosphoprotein"
  }
  
  # ------------ B-factor of the sulfur atom from methionine
  b <- factorB(pdbID)
  b <- b[which(b$chain == chain),]
  for (i in 1:nrow(Features)){
    for (j in 1:nrow(b)){
      if (Features$Met[i] == b$resno[j] + d){
        Features$Bfactor[i] <- b$b[j]
      }
    }
  }
  rm(b)
  
  # ---------- Computing the sulfur atom's depth 
  a <- dpx(pdbID, email)
  a <- a[which(a$chain == chain),]
  a$resno <- a$resno + d
  
  for (i in 1:nrow(Features)){
    for (j in 1:nrow(a)){
      if (Features$Met[i] == a$resno[j]){
        Features$dpx[i] <- a$dpx[j]
      }
    }
  }
  Features$dpx <- round(Features$dpx, 2)
  rm(a)
  
  # ------ Making predictions with random forests
  library('caret')
  library('randomForest')
  
  # Load predictive (trained) model (random forest)
  # load("rf_SBF.RData")
  # load("rf_log_mRMR_Int.Rda")
  load("rf_final.Rda")
  
  # Load predictors' names
  load("prdVars_log_SBF.RData")
  load("prdVars_log_ext_SBF.RData")

  
  # Load pre-process data (NA imputation)
  # load("preProc_mRMR.Rdata")
  load("preProc_final.Rdata")
  
  ######## To be predicted #######################################
  new_data <- Features
  ################################################################
  
  # Preprocess new data (NA imputation)
  # new_data_preproc <- predict(preProc_mRMR, new_data[,predVars_log_SBF])
  new_data_preproc <- predict(preProc_final, new_data[,predVars_log_SBF])

  # Log-transform
  new_data_preproc[, predVars_log_SBF] <- sapply(new_data_preproc[, predVars_log_SBF], log1p)

  # Compute the interaction terms
  new_data_preproc$SASA.pdb..NT_M <- new_data_preproc$SASA.pdb * new_data_preproc$NT_M
  new_data_preproc$NT_M..Met2S <- new_data_preproc$NT_M * new_data_preproc$Met2S
  new_data_preproc$NT_M..CT_F <- new_data_preproc$NT_M * new_data_preproc$CT_F
  new_data_preproc$SASA.pdb..Met2Y <- new_data_preproc$SASA.pdb * new_data_preproc$Met2Y

  new_data_preproc <- new_data_preproc[, predVars_log_ext_SBF]
  
  # Predict probabilitiy of oxidation
  # new_data_pred_prob <- predict(rf_log_mRMR_Int, new_data_preproc[,predVars_log_ext_SBF], type = "prob")[,1]
  new_data_pred_prob <- predict(rf_final, new_data_preproc[,predVars_log_ext_SBF], type = "prob")[,1]

  Features$probability <- new_data_pred_prob
  
  new_data_pred <- factor(ifelse(new_data_pred_prob > cut,"Yes", "No"),
                          levels = c('Yes', 'No'))
  Features$prediction <- new_data_pred
  
  save(Features, file=paste("", pdbID, "_MOPM.Rda", sep=""))
  
  print("WORK DONE. CHECK YOUR RESULT FILE:")
  print(paste(pdbID, "_MOPM.Rda", sep=""))
  
  output <- Features[,c(1, dim(Features)[2], dim(Features)[2]-1)]

  ## ------ Tidying up the working directory
  sub <- paste(pdbID, ".pdb", sep="")
  com <-  paste(pdbID, "_", chain, ".pdb", sep="")
  system(paste("rm", sub)) # Unix
  system(paste("rm", com)) # Unix
  system("rm -rf split_chain") # Unix
  # system("cmd.exe", input = paste("del", sub)) # Windows
  # system("cmd.exe", input = paste("del", com)) # Unix
  # system("cmd.exe", input = "rd /s /q split_chain") # Unix

  ## ---------------------------------------

  return(output)
}