# Script on the CpG mutations
# Run grep command to get the previous and following position of each candidates on the BP_resolution vcf files of each sample. 
# Formate into multiple tables as: col 1 chromosome, col 2 position, col 3 base.  

# Import samples name:
samples <- read.csv("pedigree.ped", sep =" ", header=FALSE) # From previous analysis

# Per categories:
# CtoA
CtoA_CpG=0
CtoA_no_CpG=0

for (s in samples[,3]){
    if(file.size(paste0("pos_", s, "_CtoA.txt"))==0){
    }else{
    pos <- read.csv(paste0("pos_", s, "_CtoA.txt"), sep ="\t", header=FALSE)
    data <- read.csv(paste0("CpG_all_", s, "_CtoA_used.txt"), sep ="\t", header=FALSE)
    for (i in 1:nrow(pos)){
        raw_id=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2])
        if(data[raw_id,3]=="G"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]-1)
            if(data[raw_id_n,3]=="C"){
                CtoA_CpG=CtoA_CpG+1
            }else{CtoA_no_CpG=CtoA_no_CpG+1}
        }else if (data[raw_id,3]=="C"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]+1)
            if(data[raw_id_n,3]=="G"){
                CtoA_CpG=CtoA_CpG+1
            }else{CtoA_no_CpG=CtoA_no_CpG+1}
        }
    }
}}

# CtoG
CtoG_CpG=0
CtoG_no_CpG=0

for (s in samples[,3]){
    if(file.size(paste0("pos_", s, "_CtoG.txt"))==0){}
    else{
    pos <- read.csv(paste0("pos_", s, "_CtoG.txt"), sep ="\t", header=FALSE)
    data <- read.csv(paste0("CpG_all_", s, "_CtoG_used.txt"), sep ="\t", header=FALSE)
    for (i in 1:nrow(pos)){
        raw_id=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2])
        if(data[raw_id,3]=="G"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]-1)
            if(data[raw_id_n,3]=="C"){
                CtoG_CpG=CtoG_CpG+1
            }else{CtoG_no_CpG=CtoG_no_CpG+1}
        }else if (data[raw_id,3]=="C"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]+1)
            if(data[raw_id_n,3]=="G"){
                CtoG_CpG=CtoG_CpG+1
            }else{CtoG_no_CpG=CtoG_no_CpG+1}
        }
    }
}}

# CtoT
CtoT_CpG=0
CtoT_no_CpG=0

for (s in samples[,3]){
    if(file.size(paste0("pos_", s, "_CtoT.txt"))==0){}
    else{
    pos <- read.csv(paste0("pos_", s, "_CtoT.txt"), sep ="\t", header=FALSE)
    data <- read.csv(paste0("CpG_all_", s, "_CtoT_used.txt"), sep ="\t", header=FALSE)
    for (i in 1:nrow(pos)){
        raw_id=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2])
        if(data[raw_id,3]=="G"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]-1)
            if(data[raw_id_n,3]=="C"){
                CtoT_CpG=CtoT_CpG+1
            }else{CtoT_no_CpG=CtoT_no_CpG+1}
        }else if (data[raw_id,3]=="C"){
            raw_id_n=which(as.character(data[,1])==as.character(pos[i,1]) & data[,2]==pos[i,2]+1)
            if(data[raw_id_n,3]=="G"){
                CtoT_CpG=CtoT_CpG+1
            }else{CtoT_no_CpG=CtoT_no_CpG+1}
        }
    }
}}

print("CtoA, non CpG and CpG")
CtoA_CpG+CtoA_no_CpG
CtoA_no_CpG
CtoA_CpG
print("CtoG, non CpG and CpG")
CtoG_CpG+CtoG_no_CpG
CtoG_no_CpG
CtoG_CpG
print("CtoT, non CpG and CpG")
CtoT_CpG+CtoT_no_CpG
CtoT_no_CpG
CtoT_CpG

