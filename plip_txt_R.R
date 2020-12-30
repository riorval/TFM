#función a la que se le introduce un archivo txt
sacar_interacciones <- function(pliptxt){
  #Leer las lineas del archivo
  plip1 <- readLines(pliptxt)
  #conseguir el nombre del ligando
  library(stringr)
  fila_ligando <- str_split(plip1[1], "_")
  ligando <- unlist(fila_ligando)[3]

  #encontrar filas con tipos de interacción
  interac_types <- plip1[grep("RESNR", plip1, value = FALSE)-2][-1]
  lines_interac_types <- (grep("RESNR", plip1, value = FALSE)-2)[-1]

  #encontrar filas con info del ligando
  LIG_UNL <- grep("LIG|UNL", plip1, value = FALSE)
  #quitar filas metal
  LIGnometal <- LIG_UNL[which(LIG_UNL>lines_interac_types[1])]

  interacciones <- data.frame()
  for (i in 1:length(lines_interac_types)){
    inter_type <- interac_types[i]
    if (i < length(lines_interac_types)){
      filas_interac <- LIGnometal[which(LIGnometal>lines_interac_types[i]
                                        & LIGnometal<lines_interac_types[i+1])[-1]]
    } else {
      filas_interac <- LIGnometal[which(LIGnometal>lines_interac_types[i])[-1]]
    }

    #reiniciar vector y data frame
    info <- data.frame()
    residuo <- vector()
    #para cada tipo de interacción
    for (i in 1:length(filas_interac)) {
      resnum <- unlist(str_split(plip1[filas_interac[i]], " "))[2]
      restype <- unlist(str_split(plip1[filas_interac[i]], " "))[6]
      residuo[i] <- paste0(resnum,restype)
    }
    info <- data.frame(ligando, inter_type, residuo)
    interacciones <- rbind(interacciones, info)
  }
  return(interacciones)
}

#crear data frame vacío

interaccionesTOT <- data.frame()

Batch_interac <- function(carpeta){
  #Crear una lista de archivos
  archivosplip <- list.files(carpeta, pattern = ".txt")
  #crear data frame vacío
  interacs_file <- data.frame()
  for (i in 1:length(archivosplip)){
    interacs_file <- sacar_interacciones(archivosplip[i])
    interaccionesTOT <- rbind(interaccionesTOT, interacs_file)
  }
  colnames(interaccionesTOT) <- c("Ligando", "Tipo_interac", "Residuo")
  return(interaccionesTOT)
}


#Para sacar las interacciones de 6n2w

plip1 <- readLines("6n2w_NDGA_5loxchainB_plip.txt")
interac_types <- plip1[grep("RESNR", plip1, value = FALSE)-2][-4]
lines_interac_types <- (grep("RESNR", plip1, value = FALSE)-2)[-4]
#encontrar filas con info del ligando
LIG_30Z <- grep("30Z", plip1, value = FALSE)
#repartir filas de ligando en diferentes tipos de interacción
interacciones <- data.frame()
for (i in 1:length(lines_interac_types)){
  inter_type <- interac_types[i]
  if (i < length(lines_interac_types)){
    filas_interac <- LIG_30Z[which(LIG_30Z>lines_interac_types[i]
                                      & LIG_30Z<lines_interac_types[i+1])]
  } else {
    filas_interac <- LIG_30Z[which(LIG_30Z>lines_interac_types[i])]
  }
  #reiniciar vector y data frame
  info <- data.frame()
  residuo <- vector()
  #para cada tipo de interacción
  for (i in 1:length(filas_interac)) {
    resnum <- unlist(str_split(plip1[filas_interac[i]], " "))[2]
    restype <- unlist(str_split(plip1[filas_interac[i]], " "))[6]
    residuo[i] <- paste0(resnum,restype)
  }
  ligando <- paste0("30Z_6n2w")
  info <- data.frame(ligando, inter_type, residuo)
  interacciones <- rbind(interacciones, info)
}

