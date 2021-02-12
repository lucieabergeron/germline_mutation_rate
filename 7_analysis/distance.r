# 1. Simulate random mutations ###########################################################################################################
tot_size= sum(225584828,204787373,185818997,172585720,
            190429646,180051392,169600520,144306982,
            129882849,92844088,133663169,125506784,
            108979918,127894412,111343173,77216781,
            95684472,70235451,53671032,74971481)

size1 <- 225584828
size2 <- size1 + 204787373
size3 <- size2 + 185818997
size4 <- size3 + 172585720
size5 <- size4 + 190429646
size6 <- size5 + 180051392
size7 <- size6 + 169600520
size8 <- size7 + 144306982
size9 <- size8 + 129882849
size10 <- size9 + 92844088
size11 <- size10 + 133663169
size12 <- size11 + 125506784
size13 <- size12 + 108979918
size14 <- size13 + 127894412
size15 <- size14 + 111343173
size16 <- size15 + 77216781
size17 <- size16 + 95684472
size18 <- size17 + 70235451
size19 <- size18 + 53671032
size20 <- size19 + 74971481

# 10 simulations
for(s in 1:10){
  fake_mut=runif(663, 1, tot_size)
  fake_mut_1=fake_mut[which(fake_mut>=1 & fake_mut<=size1)]
  fake_mut_2=fake_mut[which(fake_mut>size1 & fake_mut<=size2)]
  fake_mut_3=fake_mut[which(fake_mut>size2 & fake_mut<=size3)]
  fake_mut_4=fake_mut[which(fake_mut>size3 & fake_mut<=size4)]
  fake_mut_5=fake_mut[which(fake_mut>size4 & fake_mut<=size5)]
  fake_mut_6=fake_mut[which(fake_mut>size5 & fake_mut<=size6)]
  fake_mut_7=fake_mut[which(fake_mut>size6 & fake_mut<=size7)]
  fake_mut_8=fake_mut[which(fake_mut>size7 & fake_mut<=size8)]
  fake_mut_9=fake_mut[which(fake_mut>size8 & fake_mut<=size9)]
  fake_mut_10=fake_mut[which(fake_mut>size9 & fake_mut<=size10)]
  fake_mut_11=fake_mut[which(fake_mut>size10 & fake_mut<=size11)]
  fake_mut_12=fake_mut[which(fake_mut>size11 & fake_mut<=size12)]
  fake_mut_13=fake_mut[which(fake_mut>size12 & fake_mut<=size13)]
  fake_mut_14=fake_mut[which(fake_mut>size13 & fake_mut<=size14)]
  fake_mut_15=fake_mut[which(fake_mut>size14 & fake_mut<=size15)]
  fake_mut_16=fake_mut[which(fake_mut>size15 & fake_mut<=size16)]
  fake_mut_17=fake_mut[which(fake_mut>size16 & fake_mut<=size17)]
  fake_mut_18=fake_mut[which(fake_mut>size17 & fake_mut<=size18)]
  fake_mut_19=fake_mut[which(fake_mut>size18 & fake_mut<=size19)]
  fake_mut_20=fake_mut[which(fake_mut>size19 & fake_mut<=size20)]

  fake_pos_list <- c("fake_mut_1","fake_mut_2","fake_mut_3","fake_mut_4","fake_mut_5","fake_mut_6","fake_mut_7","fake_mut_8","fake_mut_9","fake_mut_10","fake_mut_11","fake_mut_12","fake_mut_13","fake_mut_14","fake_mut_15","fake_mut_16","fake_mut_17","fake_mut_18","fake_mut_19","fake_mut_20")


  for(j in fake_pos_list){
      distance=vector()
      sorted=sort(get(j))
      for(i in 2:length(get(j))){
          dist=sorted[i]-sorted[i-1]
          distance=cbind(distance,as.numeric(dist))
      }
      assign(paste0("distance_",j),distance)
  }

  distance_all_fake=c(distance_fake_mut_1,distance_fake_mut_2,distance_fake_mut_3,distance_fake_mut_4,distance_fake_mut_5,distance_fake_mut_6,distance_fake_mut_7,distance_fake_mut_8,distance_fake_mut_9,distance_fake_mut_10,distance_fake_mut_11,distance_fake_mut_12,distance_fake_mut_13,distance_fake_mut_14,distance_fake_mut_15,distance_fake_mut_16,distance_fake_mut_17,distance_fake_mut_18,distance_fake_mut_19,distance_fake_mut_20)

  # Create a table of fake mut:
  line1=cbind("chr1", as.integer(fake_mut_1))
  line2=cbind("chr2", as.integer(fake_mut_2-size1))
  line3=cbind("chr3", as.integer(fake_mut_3-size2))
  line4=cbind("chr4", as.integer(fake_mut_4-size3))
  line5=cbind("chr5", as.integer(fake_mut_5-size4))
  line6=cbind("chr6", as.integer(fake_mut_6-size5))
  line7=cbind("chr7", as.integer(fake_mut_7-size6))
  line8=cbind("chr8", as.integer(fake_mut_8-size7))
  line9=cbind("chr9", as.integer(fake_mut_9-size8))
  line10=cbind("chr10", as.integer(fake_mut_10-size9))
  line11=cbind("chr11", as.integer(fake_mut_11-size10))
  line12=cbind("chr12", as.integer(fake_mut_12-size11))
  line13=cbind("chr13", as.integer(fake_mut_13-size12))
  line14=cbind("chr14", as.integer(fake_mut_14-size13))
  line15=cbind("chr15", as.integer(fake_mut_15-size14))
  line16=cbind("chr16", as.integer(fake_mut_16-size15))
  line17=cbind("chr17", as.integer(fake_mut_17-size16))
  line18=cbind("chr18", as.integer(fake_mut_18-size17))
  line19=cbind("chr19", as.integer(fake_mut_19-size18))
  line20=cbind("chr20", as.integer(fake_mut_20-size19))

  table=rbind(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12,line13,line14,line15,line16,line17,line18,line19,line20)
  colnames(table) <- c("chrom", "pos")
  write.table(table, paste0("fake_mut_",s,".txt"),row.names=FALSE)
}

# 2. Distance between mutation and qq plot #################################################################################################################
library(plyr) # to have count()
# Import data from detection of de novo candidates
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
# direct = 'directory with the candidate mutations'

chr_size= c(225584828, 204787373,185818997,172585720,
            190429646,180051392,169600520,144306982,
            129882849,92844088,133663169,125506784,
            108979918,127894412,111343173,77216781,
            95684472,70235451,53671032,74971481)
tot=chr_size[1]+chr_size[2]+chr_size[3]+chr_size[4]+chr_size[5]+chr_size[6]+chr_size[7]+chr_size[8]+chr_size[9]+chr_size[10]+chr_size[11]+chr_size[12]+chr_size[13]+chr_size[14]+chr_size[15]+chr_size[16]+chr_size[17]+chr_size[18]+chr_size[19]+chr_size[20]
lab_name=seq(1,20,1)

distance_between_indiv = vector()
distance_all=list()
list_indiv = c("data_R16089", "data_R16053", "data_R16062", "data_R15089", "data_R15107", "data_R15123", "data_R16034", "data_R15140", "data_R14024", "data_R12006", "data_R13024", "data_R13074", "data_R13023", "data_R16100", "data_R17015", "data_R17061", "data_R17032", "data_R17037", "data_R16104")

for(indiv in list_indiv){
    data_denovo=get(indiv)[,1:8]

# Find position per chromosomes
pos_chr1=vector()
pos_chr2=vector()
pos_chr3=vector()
pos_chr4=vector()
pos_chr5=vector()
pos_chr6=vector()
pos_chr7=vector()
pos_chr8=vector()
pos_chr9=vector()
pos_chr10=vector()
pos_chr11=vector()
pos_chr12=vector()
pos_chr13=vector()
pos_chr14=vector()
pos_chr15=vector()
pos_chr16=vector()
pos_chr17=vector()
pos_chr18=vector()
pos_chr19=vector()
pos_chr20=vector()
for(i in seq(1,nrow(data_denovo))){
  if(data_denovo[i,1]=='chr1'){pos_chr1=c(pos_chr1,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr2'){pos_chr2=c(pos_chr2,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr3'){pos_chr3=c(pos_chr3,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr4'){pos_chr4=c(pos_chr4,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr5'){pos_chr5=c(pos_chr5,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr6'){pos_chr6=c(pos_chr6,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr7'){pos_chr7=c(pos_chr7,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr8'){pos_chr8=c(pos_chr8,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr9'){pos_chr9=c(pos_chr9,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr10'){pos_chr10=c(pos_chr10,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr11'){pos_chr11=c(pos_chr11,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr12'){pos_chr12=c(pos_chr12,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr13'){pos_chr13=c(pos_chr13,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr14'){pos_chr14=c(pos_chr14,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr15'){pos_chr15=c(pos_chr15,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr16'){pos_chr16=c(pos_chr16,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr17'){pos_chr17=c(pos_chr17,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr18'){pos_chr18=c(pos_chr18,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr19'){pos_chr19=c(pos_chr19,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr20'){pos_chr20=c(pos_chr20,data_denovo[i,2])
  }
}

# Find distance for each chromosomes
pos_list <- c("pos_chr1","pos_chr2","pos_chr3","pos_chr4","pos_chr5","pos_chr6","pos_chr7","pos_chr8","pos_chr9","pos_chr10","pos_chr11","pos_chr12","pos_chr13","pos_chr14","pos_chr15","pos_chr16","pos_chr17","pos_chr18","pos_chr19","pos_chr20")
rm(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr12,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)
for(j in pos_list){
    distance=vector()
    sorted=sort(get(j))
    for(i in 2:length(get(j))){
        dist=sorted[i]-sorted[i-1]
        distance=c(distance,dist)
    }
    assign(paste0("distance_",j),distance)
}

distance_all=c(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr13,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)
distance_all <- distance_all[!is.na(distance_all)]

distance_between_indiv=c(distance_between_indiv,distance_all)
}
distance_within_same_ind = distance_between_indiv

# Between related individuals 2 by 2 #################################################################
data_r1=cbind(rbind(data_R14024[,1:2], data_R12006[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R12006', nrow(data_R12006))))
data_r2=cbind(rbind(data_R14024[,1:2], data_R13024[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R13024', nrow(data_R13024))))
data_r3=cbind(rbind(data_R14024[,1:2], data_R16100[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16100', nrow(data_R16100))))
data_r4=cbind(rbind(data_R14024[,1:2], data_R13074[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R13074', nrow(data_R13074))))
data_r5=cbind(rbind(data_R14024[,1:2], data_R13023[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R13023', nrow(data_R13023))))
data_r6=cbind(rbind(data_R12006[,1:2], data_R13024[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R13024', nrow(data_R13024))))
data_r7=cbind(rbind(data_R12006[,1:2], data_R16100[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16100', nrow(data_R16100))))
data_r8=cbind(rbind(data_R12006[,1:2], data_R13074[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R13074', nrow(data_R13074))))
data_r9=cbind(rbind(data_R12006[,1:2], data_R13023[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R13023', nrow(data_R13023))))
data_r10=cbind(rbind(data_R13024[,1:2], data_R16100[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16100', nrow(data_R16100))))
data_r11=cbind(rbind(data_R13024[,1:2], data_R13074[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R13074', nrow(data_R13074))))
data_r12=cbind(rbind(data_R13024[,1:2], data_R13023[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R13023', nrow(data_R13023))))
data_r13=cbind(rbind(data_R16100[,1:2], data_R13074[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R13074', nrow(data_R13074))))
data_r14=cbind(rbind(data_R16100[,1:2], data_R13023[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R13023', nrow(data_R13023))))
data_r15=cbind(rbind(data_R13074[,1:2], data_R13023[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R13023', nrow(data_R13023))))
data_r16=cbind(rbind(data_R13074[,1:2], data_R17032[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R17032', nrow(data_R17032))))
data_r17=cbind(rbind(data_R13023[,1:2], data_R16104[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R16104', nrow(data_R16104))))
data_r18=cbind(rbind(data_R16100[,1:2], data_R17015[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R17015', nrow(data_R17015))))
data_r19=cbind(rbind(data_R16100[,1:2], data_R17061[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R17061', nrow(data_R17061))))
data_r20=cbind(rbind(data_R16100[,1:2], data_R17037[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R17037', nrow(data_R17037))))
data_r21=cbind(rbind(data_R16100[,1:2], data_R17032[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R17032', nrow(data_R17032))))
data_r22=cbind(rbind(data_R16100[,1:2], data_R16104[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R16104', nrow(data_R16104))))
data_r23=cbind(rbind(data_R17015[,1:2], data_R17061[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R17061', nrow(data_R17061))))
data_r24=cbind(rbind(data_R17015[,1:2], data_R17037[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R17037', nrow(data_R17037))))
data_r25=cbind(rbind(data_R17015[,1:2], data_R17032[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R17032', nrow(data_R17032))))
data_r26=cbind(rbind(data_R17015[,1:2], data_R16104[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R16104', nrow(data_R16104))))
data_r27=cbind(rbind(data_R17061[,1:2], data_R17037[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R17037', nrow(data_R17037))))
data_r28=cbind(rbind(data_R17061[,1:2], data_R17032[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R17032', nrow(data_R17032))))
data_r29=cbind(rbind(data_R17061[,1:2], data_R16104[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R16104', nrow(data_R16104))))
data_r30=cbind(rbind(data_R17037[,1:2], data_R17032[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R17032', nrow(data_R17032))))
data_r31=cbind(rbind(data_R17037[,1:2], data_R16104[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R16104', nrow(data_R16104))))
data_r32=cbind(rbind(data_R17032[,1:2], data_R16104[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R16104', nrow(data_R16104))))
data_r33=cbind(rbind(data_R15107[,1:2], data_R16089[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R16089', nrow(data_R16089))))
data_r34=cbind(rbind(data_R15107[,1:2], data_R16062[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R16062', nrow(data_R16062))))
data_r35=cbind(rbind(data_R15107[,1:2], data_R15123[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R15123', nrow(data_R15123))))
data_r36=cbind(rbind(data_R15107[,1:2], data_R16034[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R16034', nrow(data_R16034))))
data_r37=cbind(rbind(data_R15107[,1:2], data_R15140[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R15140', nrow(data_R15140))))
data_r38=cbind(rbind(data_R15107[,1:2], data_R16053[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R16053', nrow(data_R16053))))
data_r39=cbind(rbind(data_R15107[,1:2], data_R15089[,1:2]), c(rep('R15107', nrow(data_R15107)), rep('R15089', nrow(data_R15089))))
data_r40=cbind(rbind(data_R16089[,1:2], data_R16062[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R16062', nrow(data_R16062))))
data_r41=cbind(rbind(data_R16089[,1:2], data_R15123[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R15123', nrow(data_R15123))))
data_r42=cbind(rbind(data_R16089[,1:2], data_R16034[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R16034', nrow(data_R16034))))
data_r43=cbind(rbind(data_R16089[,1:2], data_R15140[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R15140', nrow(data_R15140))))
data_r44=cbind(rbind(data_R16089[,1:2], data_R16053[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R16053', nrow(data_R16053))))
data_r45=cbind(rbind(data_R16089[,1:2], data_R15089[,1:2]), c(rep('R16089', nrow(data_R16089)), rep('R15089', nrow(data_R15089))))
data_r46=cbind(rbind(data_R16062[,1:2], data_R15123[,1:2]), c(rep('R16062', nrow(data_R16062)), rep('R15123', nrow(data_R15123))))
data_r47=cbind(rbind(data_R16062[,1:2], data_R16034[,1:2]), c(rep('R16062', nrow(data_R16062)), rep('R16034', nrow(data_R16034))))
data_r48=cbind(rbind(data_R16062[,1:2], data_R15140[,1:2]), c(rep('R16062', nrow(data_R16062)), rep('R15140', nrow(data_R15140))))
data_r49=cbind(rbind(data_R16062[,1:2], data_R16053[,1:2]), c(rep('R16062', nrow(data_R16062)), rep('R16053', nrow(data_R16053))))
data_r50=cbind(rbind(data_R16062[,1:2], data_R15089[,1:2]), c(rep('R16062', nrow(data_R16062)), rep('R15089', nrow(data_R15089))))
data_r51=cbind(rbind(data_R15123[,1:2], data_R16034[,1:2]), c(rep('R15123', nrow(data_R15123)), rep('R16034', nrow(data_R16034))))
data_r52=cbind(rbind(data_R15123[,1:2], data_R15140[,1:2]), c(rep('R15123', nrow(data_R15123)), rep('R15140', nrow(data_R15140))))
data_r53=cbind(rbind(data_R15123[,1:2], data_R16053[,1:2]), c(rep('R15123', nrow(data_R15123)), rep('R16053', nrow(data_R16053))))
data_r54=cbind(rbind(data_R15123[,1:2], data_R15089[,1:2]), c(rep('R15123', nrow(data_R15123)), rep('R15089', nrow(data_R15089))))
data_r55=cbind(rbind(data_R16034[,1:2], data_R15140[,1:2]), c(rep('R16034', nrow(data_R16034)), rep('R15140', nrow(data_R15140))))
data_r56=cbind(rbind(data_R16034[,1:2], data_R16053[,1:2]), c(rep('R16034', nrow(data_R16034)), rep('R16053', nrow(data_R16053))))
data_r57=cbind(rbind(data_R16034[,1:2], data_R15089[,1:2]), c(rep('R16034', nrow(data_R16034)), rep('R15089', nrow(data_R15089))))
data_r58=cbind(rbind(data_R15140[,1:2], data_R16053[,1:2]), c(rep('R15140', nrow(data_R15140)), rep('R16053', nrow(data_R16053))))
data_r59=cbind(rbind(data_R15140[,1:2], data_R15089[,1:2]), c(rep('R15140', nrow(data_R15140)), rep('R15089', nrow(data_R15089))))
data_r60=cbind(rbind(data_R16053[,1:2], data_R15089[,1:2]), c(rep('R16053', nrow(data_R16053)), rep('R15089', nrow(data_R15089))))
#
distance_related_indiv = vector()
distance_all=list()
list_data = c("data_r1", "data_r2", "data_r3", "data_r4", "data_r5", "data_r6", "data_r7", "data_r8", "data_r9", "data_r10", "data_r11", "data_r12", "data_r13", "data_r14", "data_r15", "data_r16", "data_r17", "data_r18", "data_r19", "data_r20", "data_r21", "data_r22", "data_r23", "data_r24", "data_r25", "data_r26", "data_r27", "data_r28", "data_r29", "data_r30", "data_r31", "data_r32", "data_r33", "data_r34", "data_r35", "data_r36", "data_r37", "data_r38", "data_r39", "data_r40","data_r41","data_r42", "data_r43", "data_r44", "data_r45", "data_r46", "data_r47", "data_r48", "data_r49", "data_r50", "data_r51", "data_r52", "data_r53", "data_r54", "data_r55", "data_r56", "data_r57", "data_r58", "data_r59", "data_r60")

for(indiv in list_data){
    data_denovo=get(indiv)[,1:3]
# Distance
pos_chr1=vector()
pos_chr2=vector()
pos_chr3=vector()
pos_chr4=vector()
pos_chr5=vector()
pos_chr6=vector()
pos_chr7=vector()
pos_chr8=vector()
pos_chr9=vector()
pos_chr10=vector()
pos_chr11=vector()
pos_chr12=vector()
pos_chr13=vector()
pos_chr14=vector()
pos_chr15=vector()
pos_chr16=vector()
pos_chr17=vector()
pos_chr18=vector()
pos_chr19=vector()
pos_chr20=vector()
for(i in seq(1,nrow(data_denovo))){
  if(data_denovo[i,1]=='chr1'){pos_chr1=rbind(pos_chr1,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr2'){pos_chr2=rbind(pos_chr2,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr3'){pos_chr3=rbind(pos_chr3,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr4'){pos_chr4=rbind(pos_chr4,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr5'){pos_chr5=rbind(pos_chr5,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr6'){pos_chr6=rbind(pos_chr6,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr7'){pos_chr7=rbind(pos_chr7,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr8'){pos_chr8=rbind(pos_chr8,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr9'){pos_chr9=rbind(pos_chr9,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr10'){pos_chr10=rbind(pos_chr10,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr11'){pos_chr11=rbind(pos_chr11,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr12'){pos_chr12=rbind(pos_chr12,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr13'){pos_chr13=rbind(pos_chr13,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr14'){pos_chr14=rbind(pos_chr14,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr15'){pos_chr15=rbind(pos_chr15,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr16'){pos_chr16=rbind(pos_chr16,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr17'){pos_chr17=rbind(pos_chr17,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr18'){pos_chr18=rbind(pos_chr18,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr19'){pos_chr19=rbind(pos_chr19,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr20'){pos_chr20=rbind(pos_chr20,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }
}

pos_list <- c("pos_chr1","pos_chr2","pos_chr3","pos_chr4","pos_chr5","pos_chr6","pos_chr7","pos_chr8","pos_chr9","pos_chr10","pos_chr11","pos_chr12","pos_chr13","pos_chr14","pos_chr15","pos_chr16","pos_chr17","pos_chr18","pos_chr19","pos_chr20")
rm(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr12,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)
for(j in pos_list){
    distance=vector()
    set=get(j)
    if(is.null(dim(set))){
        assign(paste0("distance_",j), NA)
    }else{
        sorted=set[order(set[,1]),]
        for(i in 2:nrow(set)){
            if(nrow(sorted)<2){
              distance=c(distance,NA)
            }else{
                if(sorted[i,2]==sorted[i-1,2]){
                    distance=c(distance,NA)
                }else{
                 dist=sorted[i,1]-sorted[i-1,1]
                distance=c(distance,dist)
               }
            }
        }
    assign(paste0("distance_",j),distance)
    }
}

distance_all=c(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr13,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)
distance_all <- distance_all[!is.na(distance_all)]
distance_related_indiv=c(distance_related_indiv,distance_all)

}

# 2 by 2 between non related ##
data_nr1=cbind(rbind(data_R14024[,1:2], data_R17015[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R17015', nrow(data_R17015))))
data_nr2=cbind(rbind(data_R14024[,1:2], data_R17061[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R17061', nrow(data_R17061))))
data_nr3=cbind(rbind(data_R14024[,1:2], data_R17037[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R17037', nrow(data_R17037))))
data_nr4=cbind(rbind(data_R14024[,1:2], data_R17032[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R17032', nrow(data_R17032))))
data_nr5=cbind(rbind(data_R14024[,1:2], data_R16104[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16104', nrow(data_R16104))))
data_nr6=cbind(rbind(data_R14024[,1:2], data_R15107[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R15107', nrow(data_R15107))))
data_nr7=cbind(rbind(data_R14024[,1:2], data_R16089[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16089', nrow(data_R16089))))
data_nr8=cbind(rbind(data_R14024[,1:2], data_R16062[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16062', nrow(data_R16062))))
data_nr9=cbind(rbind(data_R14024[,1:2], data_R15123[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R15123', nrow(data_R15123))))
data_nr10=cbind(rbind(data_R14024[,1:2], data_R16034[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16034', nrow(data_R16034))))
data_nr11=cbind(rbind(data_R14024[,1:2], data_R15140[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R15140', nrow(data_R15140))))
data_nr12=cbind(rbind(data_R14024[,1:2], data_R16053[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R16053', nrow(data_R16053))))
data_nr13=cbind(rbind(data_R14024[,1:2], data_R15089[,1:2]), c(rep('R14024', nrow(data_R14024)), rep('R15089', nrow(data_R15089))))
data_nr14=cbind(rbind(data_R12006[,1:2], data_R17061[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R17061', nrow(data_R17061))))
data_nr15=cbind(rbind(data_R12006[,1:2], data_R17037[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R17037', nrow(data_R17037))))
data_nr16=cbind(rbind(data_R12006[,1:2], data_R17032[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R17032', nrow(data_R17032))))
data_nr17=cbind(rbind(data_R12006[,1:2], data_R16104[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16104', nrow(data_R16104))))
data_nr18=cbind(rbind(data_R12006[,1:2], data_R15107[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R15107', nrow(data_R15107))))
data_nr19=cbind(rbind(data_R12006[,1:2], data_R16089[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16089', nrow(data_R16089))))
data_nr20=cbind(rbind(data_R12006[,1:2], data_R16062[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16062', nrow(data_R16062))))
data_nr21=cbind(rbind(data_R12006[,1:2], data_R15123[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R15123', nrow(data_R15123))))
data_nr22=cbind(rbind(data_R12006[,1:2], data_R16034[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16034', nrow(data_R16034))))
data_nr23=cbind(rbind(data_R12006[,1:2], data_R15140[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R15140', nrow(data_R15140))))
data_nr24=cbind(rbind(data_R12006[,1:2], data_R16053[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R16053', nrow(data_R16053))))
data_nr25=cbind(rbind(data_R12006[,1:2], data_R15089[,1:2]), c(rep('R12006', nrow(data_R12006)), rep('R15089', nrow(data_R15089))))
data_nr26=cbind(rbind(data_R13024[,1:2], data_R17015[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R17015', nrow(data_R17015))))
data_nr27=cbind(rbind(data_R13024[,1:2], data_R17061[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R17061', nrow(data_R17061))))
data_nr28=cbind(rbind(data_R13024[,1:2], data_R17037[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R17037', nrow(data_R17037))))
data_nr29=cbind(rbind(data_R13024[,1:2], data_R17032[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R17032', nrow(data_R17032))))
data_nr30=cbind(rbind(data_R13024[,1:2], data_R16104[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16104', nrow(data_R16104))))
data_nr31=cbind(rbind(data_R13024[,1:2], data_R15107[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R15107', nrow(data_R15107))))
data_nr32=cbind(rbind(data_R13024[,1:2], data_R16089[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16089', nrow(data_R16089))))
data_nr33=cbind(rbind(data_R13024[,1:2], data_R16062[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16062', nrow(data_R16062))))
data_nr34=cbind(rbind(data_R13024[,1:2], data_R15123[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R15123', nrow(data_R15123))))
data_nr35=cbind(rbind(data_R13024[,1:2], data_R16034[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16034', nrow(data_R16034))))
data_nr36=cbind(rbind(data_R13024[,1:2], data_R15140[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R15140', nrow(data_R15140))))
data_nr37=cbind(rbind(data_R13024[,1:2], data_R16053[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R16053', nrow(data_R16053))))
data_nr38=cbind(rbind(data_R13024[,1:2], data_R15089[,1:2]), c(rep('R13024', nrow(data_R13024)), rep('R15089', nrow(data_R15089))))
data_nr39=cbind(rbind(data_R16100[,1:2], data_R13074[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R13074', nrow(data_R13074))))
data_nr40=cbind(rbind(data_R16100[,1:2], data_R13023[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R13023', nrow(data_R13023))))
data_nr41=cbind(rbind(data_R16100[,1:2], data_R15107[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R15107', nrow(data_R15107))))
data_nr42=cbind(rbind(data_R16100[,1:2], data_R16089[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R16089', nrow(data_R16089))))
data_nr43=cbind(rbind(data_R16100[,1:2], data_R16062[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R16062', nrow(data_R16062))))
data_nr44=cbind(rbind(data_R16100[,1:2], data_R15123[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R15123', nrow(data_R15123))))
data_nr45=cbind(rbind(data_R16100[,1:2], data_R16034[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R16034', nrow(data_R16034))))
data_nr46=cbind(rbind(data_R16100[,1:2], data_R15140[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R15140', nrow(data_R15140))))
data_nr47=cbind(rbind(data_R16100[,1:2], data_R16053[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R16053', nrow(data_R16053))))
data_nr48=cbind(rbind(data_R16100[,1:2], data_R15089[,1:2]), c(rep('R16100', nrow(data_R16100)), rep('R15089', nrow(data_R15089))))
data_nr49=cbind(rbind(data_R17061[,1:2], data_R13023[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R13023', nrow(data_R13023))))
data_nr50=cbind(rbind(data_R17061[,1:2], data_R15107[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R15107', nrow(data_R15107))))
data_nr51=cbind(rbind(data_R17061[,1:2], data_R16089[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R16089', nrow(data_R16089))))
data_nr52=cbind(rbind(data_R17061[,1:2], data_R16062[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R16062', nrow(data_R16062))))
data_nr53=cbind(rbind(data_R17061[,1:2], data_R15123[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R15123', nrow(data_R15123))))
data_nr54=cbind(rbind(data_R17061[,1:2], data_R16034[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R16034', nrow(data_R16034))))
data_nr55=cbind(rbind(data_R17061[,1:2], data_R15140[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R15140', nrow(data_R15140))))
data_nr56=cbind(rbind(data_R17061[,1:2], data_R16053[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R16053', nrow(data_R16053))))
data_nr57=cbind(rbind(data_R17061[,1:2], data_R15089[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R15089', nrow(data_R15089))))
data_nr58=cbind(rbind(data_R13074[,1:2], data_R17015[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R17015', nrow(data_R17015))))
data_nr59=cbind(rbind(data_R13074[,1:2], data_R16104[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R16104', nrow(data_R16104))))
data_nr60=cbind(rbind(data_R13074[,1:2], data_R15107[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R15107', nrow(data_R15107))))
data_nr61=cbind(rbind(data_R13074[,1:2], data_R16089[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R16089', nrow(data_R16089))))
data_nr62=cbind(rbind(data_R13074[,1:2], data_R16062[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R16062', nrow(data_R16062))))
data_nr63=cbind(rbind(data_R13074[,1:2], data_R15123[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R15123', nrow(data_R15123))))
data_nr64=cbind(rbind(data_R13074[,1:2], data_R16034[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R16034', nrow(data_R16034))))
data_nr65=cbind(rbind(data_R13074[,1:2], data_R15140[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R15140', nrow(data_R15140))))
data_nr66=cbind(rbind(data_R13074[,1:2], data_R16053[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R16053', nrow(data_R16053))))
data_nr67=cbind(rbind(data_R13074[,1:2], data_R15089[,1:2]), c(rep('R13074', nrow(data_R13074)), rep('R15089', nrow(data_R15089))))
data_nr68=cbind(rbind(data_R17032[,1:2], data_R13023[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R13023', nrow(data_R13023))))
data_nr69=cbind(rbind(data_R17032[,1:2], data_R15107[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R15107', nrow(data_R15107))))
data_nr70=cbind(rbind(data_R17032[,1:2], data_R16089[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R16089', nrow(data_R16089))))
data_nr71=cbind(rbind(data_R17032[,1:2], data_R16062[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R16062', nrow(data_R16062))))
data_nr72=cbind(rbind(data_R17032[,1:2], data_R15123[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R15123', nrow(data_R15123))))
data_nr73=cbind(rbind(data_R17032[,1:2], data_R16034[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R16034', nrow(data_R16034))))
data_nr74=cbind(rbind(data_R17032[,1:2], data_R15140[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R15140', nrow(data_R15140))))
data_nr75=cbind(rbind(data_R17032[,1:2], data_R16053[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R16053', nrow(data_R16053))))
data_nr76=cbind(rbind(data_R17032[,1:2], data_R15089[,1:2]), c(rep('R17032', nrow(data_R17032)), rep('R15089', nrow(data_R15089))))
data_nr77=cbind(rbind(data_R13023[,1:2], data_R17015[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R17015', nrow(data_R17015))))
data_nr78=cbind(rbind(data_R13023[,1:2], data_R17037[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R17037', nrow(data_R17037))))
data_nr79=cbind(rbind(data_R13023[,1:2], data_R15107[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R15107', nrow(data_R15107))))
data_nr80=cbind(rbind(data_R13023[,1:2], data_R16089[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R16089', nrow(data_R16089))))
data_nr81=cbind(rbind(data_R13023[,1:2], data_R16062[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R16062', nrow(data_R16062))))
data_nr82=cbind(rbind(data_R13023[,1:2], data_R15123[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R15123', nrow(data_R15123))))
data_nr83=cbind(rbind(data_R13023[,1:2], data_R16034[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R16034', nrow(data_R16034))))
data_nr84=cbind(rbind(data_R13023[,1:2], data_R15140[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R15140', nrow(data_R15140))))
data_nr85=cbind(rbind(data_R13023[,1:2], data_R16053[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R16053', nrow(data_R16053))))
data_nr86=cbind(rbind(data_R13023[,1:2], data_R15089[,1:2]), c(rep('R13023', nrow(data_R13023)), rep('R15089', nrow(data_R15089))))
data_nr87=cbind(rbind(data_R16104[,1:2], data_R15107[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R15107', nrow(data_R15107))))
data_nr88=cbind(rbind(data_R16104[,1:2], data_R16089[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R16089', nrow(data_R16089))))
data_nr89=cbind(rbind(data_R16104[,1:2], data_R16062[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R16062', nrow(data_R16062))))
data_nr90=cbind(rbind(data_R16104[,1:2], data_R15123[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R15123', nrow(data_R15123))))
data_nr91=cbind(rbind(data_R16104[,1:2], data_R16034[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R16034', nrow(data_R16034))))
data_nr92=cbind(rbind(data_R16104[,1:2], data_R15140[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R15140', nrow(data_R15140))))
data_nr93=cbind(rbind(data_R16104[,1:2], data_R16053[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R16053', nrow(data_R16053))))
data_nr94=cbind(rbind(data_R16104[,1:2], data_R15089[,1:2]), c(rep('R16104', nrow(data_R16104)), rep('R15089', nrow(data_R15089))))
data_nr95=cbind(rbind(data_R17015[,1:2], data_R15107[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R15107', nrow(data_R15107))))
data_nr96=cbind(rbind(data_R17015[,1:2], data_R16089[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R16089', nrow(data_R16089))))
data_nr97=cbind(rbind(data_R17015[,1:2], data_R16062[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R16062', nrow(data_R16062))))
data_nr98=cbind(rbind(data_R17015[,1:2], data_R15123[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R15123', nrow(data_R15123))))
data_nr99=cbind(rbind(data_R17015[,1:2], data_R16034[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R16034', nrow(data_R16034))))
data_nr100=cbind(rbind(data_R17015[,1:2], data_R15140[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R15140', nrow(data_R15140))))
data_nr101=cbind(rbind(data_R17015[,1:2], data_R16053[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R16053', nrow(data_R16053))))
data_nr102=cbind(rbind(data_R17015[,1:2], data_R15089[,1:2]), c(rep('R17015', nrow(data_R17015)), rep('R15089', nrow(data_R15089))))
data_nr103=cbind(rbind(data_R17037[,1:2], data_R15107[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R15107', nrow(data_R15107))))
data_nr104=cbind(rbind(data_R17037[,1:2], data_R16089[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R16089', nrow(data_R16089))))
data_nr105=cbind(rbind(data_R17037[,1:2], data_R16062[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R16062', nrow(data_R16062))))
data_nr106=cbind(rbind(data_R17037[,1:2], data_R15123[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R15123', nrow(data_R15123))))
data_nr107=cbind(rbind(data_R17037[,1:2], data_R16034[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R16034', nrow(data_R16034))))
data_nr108=cbind(rbind(data_R17037[,1:2], data_R15140[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R15140', nrow(data_R15140))))
data_nr109=cbind(rbind(data_R17037[,1:2], data_R16053[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R16053', nrow(data_R16053))))
data_nr110=cbind(rbind(data_R17037[,1:2], data_R15089[,1:2]), c(rep('R17037', nrow(data_R17037)), rep('R15089', nrow(data_R15089))))
data_nr111=cbind(rbind(data_R17061[,1:2], data_R13074[,1:2]), c(rep('R17061', nrow(data_R17061)), rep('R13074', nrow(data_R13074))))
#
distance_non_related_indiv = vector()
distance_all=list()
list_data = c("data_nr1", "data_nr2", "data_nr3", "data_nr4", "data_nr5", "data_nr6", "data_nr7", "data_nr8", "data_nr9", "data_nr10", "data_nr11", "data_nr12", "data_nr13", "data_nr14", "data_nr15", "data_nr16", "data_nr17", "data_nr18", "data_nr19", "data_nr20", "data_nr21", "data_nr22", "data_nr23", "data_nr24", "data_nr25", "data_nr26", "data_nr27", "data_nr28", "data_nr29", "data_nr30", "data_nr31", "data_nr32", "data_nr33", "data_nr34", "data_nr35", "data_nr36", "data_nr37", "data_nr38", "data_nr39", "data_nr40","data_nr41","data_nr42", "data_nr43", "data_nr44", "data_nr45", "data_nr46", "data_nr47", "data_nr48", "data_nr49", "data_nr50", "data_nr51", "data_nr52", "data_nr53", "data_nr54", "data_nr55", "data_nr56", "data_nr57", "data_nr58", "data_nr59", "data_nr60", "data_nr70", "data_nr71", "data_nr72", "data_nr73", "data_nr74", "data_nr75", "data_nr76", "data_nr77", "data_nr78", "data_nr79", "data_nr80", "data_nr81", "data_nr82", "data_nr83", "data_nr84", "data_nr85", "data_nr86", "data_nr87", "data_nr88", "data_nr89", "data_nr90", "data_nr91", "data_nr92", "data_nr93", "data_nr94", "data_nr95", "data_nr96", "data_nr97", "data_nr98", "data_nr99", "data_nr100", "data_nr101", "data_nr102", "data_nr103", "data_nr104", "data_nr105", "data_nr106", "data_nr107", "data_nr108", "data_nr109", "data_nr110",  "data_nr111")

for(indiv in list_data){
    data_denovo=get(indiv)[,1:3]

## Distance
pos_chr1=vector()
pos_chr2=vector()
pos_chr3=vector()
pos_chr4=vector()
pos_chr5=vector()
pos_chr6=vector()
pos_chr7=vector()
pos_chr8=vector()
pos_chr9=vector()
pos_chr10=vector()
pos_chr11=vector()
pos_chr12=vector()
pos_chr13=vector()
pos_chr14=vector()
pos_chr15=vector()
pos_chr16=vector()
pos_chr17=vector()
pos_chr18=vector()
pos_chr19=vector()
pos_chr20=vector()
for(i in seq(1,nrow(data_denovo))){
  if(data_denovo[i,1]=='chr1'){pos_chr1=rbind(pos_chr1,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr2'){pos_chr2=rbind(pos_chr2,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr3'){pos_chr3=rbind(pos_chr3,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr4'){pos_chr4=rbind(pos_chr4,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr5'){pos_chr5=rbind(pos_chr5,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr6'){pos_chr6=rbind(pos_chr6,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr7'){pos_chr7=rbind(pos_chr7,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr8'){pos_chr8=rbind(pos_chr8,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr9'){pos_chr9=rbind(pos_chr9,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr10'){pos_chr10=rbind(pos_chr10,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr11'){pos_chr11=rbind(pos_chr11,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr12'){pos_chr12=rbind(pos_chr12,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr13'){pos_chr13=rbind(pos_chr13,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr14'){pos_chr14=rbind(pos_chr14,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr15'){pos_chr15=rbind(pos_chr15,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr16'){pos_chr16=rbind(pos_chr16,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr17'){pos_chr17=rbind(pos_chr17,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr18'){pos_chr18=rbind(pos_chr18,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr19'){pos_chr19=rbind(pos_chr19,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }else if(data_denovo[i,1]=='chr20'){pos_chr20=rbind(pos_chr20,data.frame(data_denovo[i,2], data_denovo[i,3]))
  }
}

pos_list <- c("pos_chr1","pos_chr2","pos_chr3","pos_chr4","pos_chr5","pos_chr6","pos_chr7","pos_chr8","pos_chr9","pos_chr10","pos_chr11","pos_chr12","pos_chr13","pos_chr14","pos_chr15","pos_chr16","pos_chr17","pos_chr18","pos_chr19","pos_chr20")
rm(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr12,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)

for(j in pos_list){
    distance=vector()
    set=get(j)
    if(is.null(dim(set))){
        assign(paste0("distance_",j), NA)
    }else{
        sorted=set[order(set[,1]),]
        for(i in 2:nrow(set)){
            if(nrow(sorted)<2){
              distance=c(distance,NA)
            }else{
                if(sorted[i,2]==sorted[i-1,2]){
                    distance=c(distance,NA)
                }else{
                 dist=sorted[i,1]-sorted[i-1,1]
                distance=c(distance,dist)
               }
            }
        }
    assign(paste0("distance_",j),distance)
    }
}

distance_all=c(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr13,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)
distance_all <- distance_all[!is.na(distance_all)]
distance_non_related_indiv=c(distance_non_related_indiv,distance_all)

}

###
## Now import the simulated for QQ_plot:
fake_mut_1 = read.csv("fake_mut_1.txt", sep=' ')

# Distance for fake_mut_1
data_denovo = fake_mut_1
pos_chr1=vector()
pos_chr2=vector()
pos_chr3=vector()
pos_chr4=vector()
pos_chr5=vector()
pos_chr6=vector()
pos_chr7=vector()
pos_chr8=vector()
pos_chr9=vector()
pos_chr10=vector()
pos_chr11=vector()
pos_chr12=vector()
pos_chr13=vector()
pos_chr14=vector()
pos_chr15=vector()
pos_chr16=vector()
pos_chr17=vector()
pos_chr18=vector()
pos_chr19=vector()
pos_chr20=vector()
for(i in seq(1,nrow(data_denovo))){
  if(data_denovo[i,1]=='chr1'){pos_chr1=c(pos_chr1,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr2'){pos_chr2=c(pos_chr2,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr3'){pos_chr3=c(pos_chr3,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr4'){pos_chr4=c(pos_chr4,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr5'){pos_chr5=c(pos_chr5,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr6'){pos_chr6=c(pos_chr6,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr7'){pos_chr7=c(pos_chr7,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr8'){pos_chr8=c(pos_chr8,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr9'){pos_chr9=c(pos_chr9,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr10'){pos_chr10=c(pos_chr10,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr11'){pos_chr11=c(pos_chr11,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr12'){pos_chr12=c(pos_chr12,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr13'){pos_chr13=c(pos_chr13,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr14'){pos_chr14=c(pos_chr14,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr15'){pos_chr15=c(pos_chr15,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr16'){pos_chr16=c(pos_chr16,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr17'){pos_chr17=c(pos_chr17,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr18'){pos_chr18=c(pos_chr18,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr19'){pos_chr19=c(pos_chr19,data_denovo[i,2])
  }else if(data_denovo[i,1]=='chr20'){pos_chr20=c(pos_chr20,data_denovo[i,2])
  }
}

## Distances:
pos_list <- c("pos_chr1","pos_chr2","pos_chr3","pos_chr4","pos_chr5","pos_chr6","pos_chr7","pos_chr8","pos_chr9","pos_chr10","pos_chr11","pos_chr12","pos_chr13","pos_chr14","pos_chr15","pos_chr16","pos_chr17","pos_chr18","pos_chr19","pos_chr20")
rm(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr13,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)

for(j in pos_list){
    distance=vector()
    sorted=sort(get(j))
    for(i in 2:length(get(j))){
        dist=sorted[i]-sorted[i-1]
        distance=c(distance,dist)
    }
    assign(paste0("distance_",j),distance)
}
distance_all=c(distance_pos_chr1,distance_pos_chr2,distance_pos_chr3,distance_pos_chr4,distance_pos_chr5,distance_pos_chr6,distance_pos_chr7,distance_pos_chr8,distance_pos_chr9,distance_pos_chr10,distance_pos_chr11,distance_pos_chr12,distance_pos_chr13,distance_pos_chr14,distance_pos_chr15,distance_pos_chr16,distance_pos_chr17,distance_pos_chr18,distance_pos_chr19,distance_pos_chr20)

distance_fake_1 = distance_all

##############################################################
# QQ plot with all:
# distance_within_same_ind, distance_related_indiv, distance_non_related_indiv
# "#9553CB", "#49C9A5", "#FF7B29"

png("QQ_plot_diff_better.png", width = 550, height = 850)
par(mar=c(10,11,4,2), mgp=c(3, 2, 0))
qqplot(log10(distance_fake_1),log10(distance_within_same_ind), xlim=c(0,9), ylim=c(0,9), pch='-', lwd=10, xlab="",ylab="",cex=3,cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n", col="#9553CB")
par(new=TRUE)
qqplot(log10(distance_fake_1),log10(distance_related_indiv), xlim=c(0,9), ylim=c(0,9), pch='-', lwd=10, xlab="",ylab="",cex=3,cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n", col="#49C9A5")
par(new=TRUE)
qqplot(log10(distance_fake_1),log10(distance_non_related_indiv), xlim=c(0,9), ylim=c(0,9), pch='-', lwd=10, xlab="",ylab="",cex=3,cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n", col="#FF7B29")
axis(1, cex.axis=3)
axis(2, cex.axis=3, las=2)
mtext(expression(paste("Expected distance (log"[10],")")), 1, line=6, cex=3)
mtext(expression(paste("Observed distance (log"[10],")")), 2, line=5, cex=3)
legend("topleft", bty = "n", fill=c("#9553CB", "#49C9A5", "#FF7B29"),c("same", "related", "non-related"), x.intersp =0.5, cex=2.8)
dev.off()

# 3. Shared between siblings:
## data table from de novo detection
data_denovo_last <- c(rep("R16089", nrow(data_R16089)), rep("R16053", nrow(data_R16053)),
                      rep("R16062", nrow(data_R16062)), rep("R15089", nrow(data_R15089)),
                      rep("R15107", nrow(data_R15107)), rep("R15123", nrow(data_R15123)),
                      rep("R16034", nrow(data_R16034)), rep("R15140", nrow(data_R15140)),
                      rep("R14024", nrow(data_R14024)), rep("R12006", nrow(data_R12006)),
                      rep("R13024", nrow(data_R13024)), rep("R13074", nrow(data_R13074)),
                      rep("R13023", nrow(data_R13023)), rep("R16100", nrow(data_R16100)),
                      rep("R17015", nrow(data_R17015)), rep("R17061", nrow(data_R17061)),
                      rep("R17032", nrow(data_R17032)), rep("R17037", nrow(data_R17037)),
                      rep("R16104", nrow(data_R16104)))


data_denovo_full <- cbind(data_denovo, data_denovo_last)
shared <- data_denovo_full[which(duplicated(data_denovo) | duplicated(data_denovo[nrow(data_denovo):1, ])[nrow(data_denovo):1]=="TRUE"),]
shared[order(shared$CHROM,shared$POS),]
write.table(shared[order(shared$CHROM,shared$POS),], "shared_denovo.tab", sep="\t")

# Chance to shared on the random data in distance
for(i in seq(1,10)){
    fake <-read.csv(paste0("../distance/fake_mut_",i,".txt"), sep =" ", header = TRUE)
    shared_fake <- fake[which(duplicated(fake) | duplicated(fake[nrow(fake):1, ])[nrow(fake):1]=="TRUE"),]
    print(shared_fake)
}


