#recoger interacciones de NDGA cristalizado

interac_NDGA6n2w <- unique(interacciones_NDGA_6n2w$residuo)

#cambiar en interacciones_inhibidores_5lox NDGA para especificar el modelo (1,2 y 3)
interacciones_inhibidores_5lox$Ligando[39:47] <- rep("NDGA_1", 9)
interacciones_inhibidores_5lox$Ligando[48:58] <- rep("NDGA_2", 11)
interacciones_inhibidores_5lox$Ligando[59:70] <- rep("NDGA_3", 12)

#de los inhibidores, separar docking NDGA del resto de inhibidores
interac_docking_NDGA <- interacciones_inhibidores_5lox[39:70,]
interac_inhibidores_noNDGA <- interacciones_inhibidores_5lox[-(39:70),]

#ver qué interacciones tienen en común NDGA dockeado con el NDGA cristalizado
interac_conservadas_NDGA <- interac_docking_NDGA[unlist(lapply(interac_NDGA6n2w,
                                                               function(x) grep(x, interac_docking_NDGA$Residuo, fixed = TRUE))),]
#se conservan las interacciones con 6/7 residuos.
#El único con el que no se conservan las interacciones en el docking de NDGA es 569PRO.

#quitar interacciones repetidas de la tabla de NDGA cristalizado
interaNoRep_NDGA_6n2w <- interacciones_NDGA_6n2w[-3,]
interaNoRep_NDGA_6n2w <- interaNoRep_NDGA_6n2w[-7,]

#comprobar que el tipo de interacción también coincide
inter_type_conservadas_NDGA <- data.frame()
for (i in 1:nrow(interaNoRep_NDGA_6n2w)){
  conservadas <- subset(interac_docking_NDGA,
                        interac_docking_NDGA$Tipo_interac == interaNoRep_NDGA_6n2w$inter_type[i]
                        & interac_docking_NDGA$Residuo == interaNoRep_NDGA_6n2w$residuo[i])
  inter_type_conservadas_NDGA <- rbind(inter_type_conservadas_NDGA, conservadas)
}

table(inter_type_conservadas_NDGA$Residuo)

#sólo 6 interacciones se conservan en los modelos de docking de NDGA
conser_norep_NDGA <- interaNoRep_NDGA_6n2w[-5,]

#(ver qué inhibidores tienen más interacciones en común con el NDGA cristalizado)
table(interac_inhib_enNDGA6n2w$Ligando)

#De las interacciones conservadas con el docking, cuáles se conservan en el resto de inhibidores conocidos
inhib_conoc_conservadas <- data.frame()
for (i in 1:nrow(conser_norep_NDGA)){
  conservadas <- subset(interac_inhibidores_noNDGA,
                        interac_inhibidores_noNDGA$Tipo_interac == conser_norep_NDGA$inter_type[i]
                        & interac_inhibidores_noNDGA$Residuo == conser_norep_NDGA$residuo[i])
  inhib_conoc_conservadas <- rbind(inhib_conoc_conservadas, conservadas)
}

table(inhib_conoc_conservadas$Residuo)
#se conservan 5 de las interacciones
#para contar en cuántos inhibidores se conserva cada interacción
inhib_conoc_conservadas_norep <- unique(inhib_conoc_conservadas)
table(inhib_conoc_conservadas_norep$Residuo)

#comprobar si la interacción no conservada en el docking de NDGA está en los inhibidores conocidos
inhibidores_569pro <- data.frame()
inhibidores_569pro <- subset(interac_inhibidores_noNDGA, interac_inhibidores_noNDGA$Residuo == "569PRO")

#crear lista de las 5 interacciones más conservadas
lista_conservadas <- conser_norep_NDGA[-4,-1]

#interacciones conservadas en 25 flavonoles
#voy a utilizar las 5 que se conservan en los inhibidores
flavon_conservadas <- data.frame()
for (i in 1:nrow(lista_conservadas)){
  conservadas <- subset(interacciones_25flavonoles_5lox,
                        interacciones_25flavonoles_5lox$Tipo_interac == lista_conservadas$inter_type[i]
                        & interacciones_25flavonoles_5lox$Residuo == lista_conservadas$residuo[i])
  flavon_conservadas <- rbind(flavon_conservadas, conservadas)
}

#elimino las interacciones repetidas
flavon_conservadas_norep <- unique(flavon_conservadas[,-1])

#cuenta de flavonoles con interacciones no repetidas
cuenta_inter_conser_flavon <- as.data.frame(table(flavon_conservadas_norep$Ligando))
#seleccionar flavonoles con más de 4 repeticiones
flavon6_4interac <- cuenta_inter_conser_flavon[which(cuenta_inter_conser_flavon$Freq>=4),]
flavon6_4interac_lista <- as.vector(flavon6_4interac[,1])

#recoger interacciones de flavonoles con más interacciones
interac_6flavonolesTODAS <- data.frame()
for (i in 1:length(flavon6_4interac_lista)){
  conservadas <- subset(interacciones_25flavonoles_5lox,
                        interacciones_25flavonoles_5lox$Ligando == flavon6_4interac_lista[i])
  interac_6flavonolesTODAS <- rbind(interac_6flavonolesTODAS, conservadas)
}
#estas son sólo las interacciones conservadas en los 6 flavonoles con más interac conservadas
interConser_6flavonoles <- data.frame()
for (i in 1:length(flavon6_4interac_lista)){
  conservadas <- subset(flavon_conservadas,
                        flavon_conservadas$Ligando == flavon6_4interac_lista[i])
  interConser_6flavonoles <- rbind(interConser_6flavonoles, conservadas)
}

#flavonoles NATURALES (9 en total)
#interacciones conservadas en flavonoles naturales
flavoNatur_conservadas <- data.frame()
for (i in 1:nrow(lista_conservadas)){
  conservadas <- subset(Interac_flavoNatur,
                        Interac_flavoNatur$Tipo_interac == lista_conservadas$inter_type[i]
                        & Interac_flavoNatur$Residuo == lista_conservadas$residuo[i])
  flavoNatur_conservadas <- rbind(flavoNatur_conservadas, conservadas)
}

#elimino las interacciones repetidas
flavoNatur_conservadas_norep <- unique(flavoNatur_conservadas)

#cuenta de flavonoles naturales con interacciones no repetidas
cuenta_inter_conser_flavoNatur <- as.data.frame(table(flavoNatur_conservadas_norep$Ligando))
#seleccionar flavonoles naturales con más de 4 repeticiones
flavoNatur_4interac <- cuenta_inter_conser_flavoNatur[which(cuenta_inter_conser_flavoNatur$Freq>=4),]
flavoNatur_4interac_lista <- as.vector(flavoNatur_4interac[,1])

#recoger interacciones de flavonoles con más interacciones
interac_3flavoNaturTODAS <- data.frame()
for (i in 1:length(flavoNatur_4interac_lista)){
  conservadas <- subset(Interac_flavoNatur,
                        Interac_flavoNatur$Ligando == flavoNatur_4interac_lista[i])
  interac_3flavoNaturTODAS <- rbind(interac_3flavoNaturTODAS, conservadas)
}
#estas son sólo las interacciones conservadas en los 6 flavonoles con más interac conservadas
interConser_3flavoNATUR <- data.frame()
for (i in 1:length(flavoNatur_4interac_lista)){
  conservadas <- subset(flavoNatur_conservadas,
                        flavoNatur_conservadas$Ligando == flavoNatur_4interac_lista[i])
  interConser_3flavoNATUR <- rbind(interConser_3flavoNATUR, conservadas)
}
