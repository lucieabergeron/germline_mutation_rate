# Script on mutation spectrum
# data obtained during de novo detection:
data_R16089 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16089.tab", sep ="\t")
data_R16053 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16053.tab", sep ="\t")
data_R16062 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16062.tab", sep ="\t")
data_R15089 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R15089.tab", sep ="\t")
data_R15107 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R15107.tab", sep ="\t")
data_R15123 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R15123.tab", sep ="\t")
data_R16034 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16034.tab", sep ="\t")
data_R15140 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R15140.tab", sep ="\t")
data_R14024 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R14024.tab", sep ="\t")
data_R12006 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R12006.tab", sep ="\t")
data_R13024 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R13024.tab", sep ="\t")
data_R13074 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R13074.tab", sep ="\t")
data_R13023 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R13023.tab", sep ="\t")
data_R16100 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16100.tab", sep ="\t")
data_R17015 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R17015.tab", sep ="\t")
data_R17061 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R17061.tab", sep ="\t")
data_R17032 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R17032.tab", sep ="\t")
data_R17037 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R17037.tab", sep ="\t")
data_R16104 = read.csv("/home/lucie/MammalianMutation/faststorage/Macaca_mulatta_33/de_novo_mutation/data_denovo_R16104.tab", sep ="\t")
# Sample:
name="all"

# Data:
data_denovo <- rbind(data_R16089[,1:8], data_R16053[,1:8], data_R16062[,1:8], data_R15089[,1:8], data_R15107[,1:8], data_R15123[,1:8], data_R16034[,1:8], data_R15140[,1:8], data_R14024[,1:8], data_R12006[,1:8], data_R13024[,1:8], data_R13074[,1:8], data_R13023[,1:8], data_R16100[,1:8], data_R17015[,1:8], data_R17061[,1:8], data_R17032[,1:8], data_R17037[,1:8], data_R16104[,1:8])
write.table(data_denovo, "all_denovo.tab", sep="\t")

# Position and name of all samples:
name=c(rep("R16089",nrow(data_R16089)),
       rep("R16053",nrow(data_R16053)),
       rep("R16062",nrow(data_R16062)),
       rep("R15089",nrow(data_R15089)),
       rep("R15107",nrow(data_R15107)),
       rep("R15123",nrow(data_R15123)),
       rep("R16034",nrow(data_R16034)),
       rep("R15140",nrow(data_R15140)),
       rep("R14024",nrow(data_R14024)),
       rep("R12006",nrow(data_R12006)),
       rep("R13024",nrow(data_R13024)),
       rep("R13074",nrow(data_R13074)),
       rep("R13023",nrow(data_R13023)),
       rep("R16100",nrow(data_R16100)),
       rep("R17015",nrow(data_R17015)),
       rep("R17061",nrow(data_R17061)),
       rep("R17032",nrow(data_R17032)),
       rep("R17037",nrow(data_R17037)),
       rep("R16104",nrow(data_R16104)))
write.table(cbind(data_denovo[,1:2],name), "all_position.tab", sep="\t")

# Caracterisation type of mutation #####################################################

AtoC=0
AtoG=0
AtoT=0
CtoG=0
CtoT=0
CtoA=0
TtoG=0
TtoC=0
TtoA=0
GtoC=0
GtoA=0
GtoT=0
for(i in seq(1,nrow(data_denovo))){
  if(data_denovo[i,6]=='A'){
    if(data_denovo[i,7]=='C'){AtoC=AtoC+1
    }else if(data_denovo[i,7]=='G'){AtoG=AtoG+1
    }else if(data_denovo[i,7]=='T'){AtoT=AtoT+1}
  }else if(data_denovo[i,6]=='C'){
    if(data_denovo[i,7]=='G'){CtoG=CtoG+1
    }else if(data_denovo[i,7]=='T'){CtoT=CtoT+1
    }else if(data_denovo[i,7]=='A'){CtoA=CtoA+1}
  }else if(data_denovo[i,6]=='T'){
    if(data_denovo[i,7]=='G'){TtoG=TtoG+1
    }else if(data_denovo[i,7]=='C'){TtoC=TtoC+1
    }else if(data_denovo[i,7]=='A'){TtoA=TtoA+1}
  }else if(data_denovo[i,6]=='G'){
    if(data_denovo[i,7]=='C'){GtoC=GtoC+1
    }else if(data_denovo[i,7]=='A'){GtoA=GtoA+1
    }else if(data_denovo[i,7]=='T'){GtoT=GtoT+1}
  }
}

#A>C (transversion) = T>G (transversion)
#A>G (transition) = T>C (transition)
#A>T (transversion) = T>A (transversion)
#C>A (transversion) = G>T (transversion)
#C>G (transversion) = G>C (transversion)
#C>T (transition) = G>A (transition)

A.C=AtoC+TtoG #T.G=TtoG+AtoC
A.G=AtoG+TtoC #T.C=TtoC+AtoG
A.T=AtoT+TtoA #T.A=TtoA+AtoT
C.A=CtoA+GtoT #G.T=GtoT+CtoA
C.G=CtoG+GtoC #G.C=GtoC+CtoG
C.T=CtoT+GtoA #G.A=GtoA+CtoT

# Information from CpG code
C.A_1=62
C.A_2=8
C.G_1=44
C.G_2=2
C.T_1=188
C.T_2=144

# Plot
png(paste0("type_of_mutation_", name, ".png"), width = 900, height = 850)
par(mar=c(10,11,4,2), mgp=c(3, 2.5, 0))
mut_type=(matrix(c(A.C,0,A.G,0,A.T,0,C.A_1,C.A_2,C.G_1,C.G_2,C.T_1,C.T_2), nrow=2))
barplot(mut_type, names.arg=c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T"), ylim=c(0,355),
        col=c("grey33", "grey72"), cex.axis = 3, cex.names=3, yaxt="n")
mtext(c("(T>G)", "(T>C)", "(T>A)", "(G>T)", "(G>C)", "(G>A)"), 1, at=c(0.7,1.9,3.1,4.3,5.5,6.7), line=5, cex=2)
#barplot(c(0,0,0,0,0,mut_type[,6]), add=TRUE, col=c("black"), density =15, cex.axis = 1.5, cex.names=1.5)
legend("topleft", bty = "n", fill=c("grey72", "gray33"),c("CpG", "non-CpG"), cex=3)
#axis(1, at=c(1,2,3,4,5,6), label=c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T"), cex.axis=3.5)
axis(2, cex.axis=3, las=2)
axis(1, at=c(0.7,1.9,3.1,4.3,5.5,6.7), label=c("", "", "", "", "", ""), cex.axis=3.5)
mtext("Type of mutation", 1, line=8, cex=3)
mtext("Number of mutations", 2, line=8, cex=3)
box()
dev.off()

#######################################################################
# Allelic balance of de novo candidates
data_ab <- c(as.character(data_R16089[,10]), as.character(data_R16053[,10]), as.character(data_R16062[,10]), as.character(data_R15089[,10]), as.character(data_R15107[,10]), as.character(data_R15123[,10]), as.character(data_R16034[,10]), as.character(data_R15140[,10]), as.character(data_R14024[,10]), as.character(data_R12006[,10]), as.character(data_R13024[,10]), as.character(data_R13074[,10]), as.character(data_R13023[,10]), as.character(data_R16100[,10]), as.character(data_R17015[,10]), as.character(data_R17061[,10]), as.character(data_R17032[,10]), as.character(data_R17037[,10]), as.character(data_R16104[,10]))

ab_all=vector()
for(i in seq(1:624)){
   alt=as.numeric(strsplit(data_ab[i],',')[[1]][2])
   tot=as.numeric(strsplit(data_ab[i],',')[[1]][1])+as.numeric(strsplit(data_ab[i],',')[[1]][2])
   ab=alt/tot
   ab_all=c(ab_all, ab)
}

png("AB_denovo.png", width = 1300, height = 850)
par(mar=c(10,16,4,2), mgp=c(3, 2, 0))
hist(ab_all, breaks=15, xlab="", ylab="", xlim=c(0,1), cex=3, cex.lab=2.5, cex.axis=2.5, main="", xaxt="n", yaxt="n")
axis(2, cex.axis=2.5, las=2)
axis(1, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), label=c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"), cex.axis=2.5)
mtext(side = 1, text = "Allelic balance" , line = 6, cex=2.5)
mtext("Frequency", 2, line=12, cex=2.5)
box()
dev.off()
