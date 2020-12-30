

smiles <- 0
ID <- 0
txtname <- 0

for (i in 1:nrow(FlavonolesChEMBL)){
  smiles[i] <- as.character(FlavonolesChEMBL[i,3])
  ID[i] <- as.character(FlavonolesChEMBL[i,1])
  txtname[i] <- paste0(ID[i], ".smi")
  write(paste(smiles[i], "  ", ID[i]), file = txtname[i])}

