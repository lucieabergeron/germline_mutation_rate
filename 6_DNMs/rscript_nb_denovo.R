# Import data ##########################################################################
data <- read.csv(paste0(direct_handling, "trio_", name, "_vcf_table.out"), sep ="\t")
########################################################################################
# We want: Offspring, Father and Mother.
# But now it is alphabetic order.
data_commun <- data[,1:6] # this have no samples name
data1 <- data[,7:12] # first sample
data2 <- data[,13:18] # second sample
data3 <- data[,19:24] # last sample
######################
# Import the pedigree information:
samples <- read.csv(paste0(direct,"vcf_tables_directories.txt"), sep =" ", header=FALSE)
samples_line=which(samples[,1]==name)
samples_off=toString(samples[samples_line,1])
samples_fa=toString(samples[samples_line,2])
samples_mo=toString(samples[samples_line,3])
samples_trio=c(samples_off, samples_fa, samples_mo)
samples_vcf=sort(samples_trio)
# Finds the correct order
order=c(which(samples_vcf==samples_off),
        which(samples_vcf==samples_fa),
        which(samples_vcf==samples_mo))
# And order the data file
data_order <- cbind(data_commun,
                    eval(parse(text=paste0("data", order[1]))),
                    eval(parse(text=paste0("data", order[2]))),
                    eval(parse(text=paste0("data", order[3]))))

########################################################################################
# Keep only the filtered ones
data_pass_site_filter <- data_order[which(data_order$FILTER=="PASS"),]
########################################################################################

# Find de novo (more strict than MV) ###################################################
denovo <- numeric(length = nrow(data_pass_site_filter))
type_denovo <- numeric(length = nrow(data_pass_site_filter))
for(j in seq(1, nrow(data_pass_site_filter))){
  # Set up what is HomRef, HomAlt, HetAlt for each line
  HomRef=paste(data_pass_site_filter[j,4], "/", data_pass_site_filter[j,4], sep="")
  alt_nucl <- strsplit(as.character(data_pass_site_filter[j,5]), split=",")[[1]]
  if(length(alt_nucl)==1){ # We have only one alt allele
    # Set homozygote Alt
    HomAlt=paste(data_pass_site_filter[j,5], "/", data_pass_site_filter[j,5], sep="")
    # Set heterozygotes
    HetAlt1=paste(data_pass_site_filter[j,4], "/", data_pass_site_filter[j,5], sep="")
    HetAlt2=paste(data_pass_site_filter[j,5], "/", data_pass_site_filter[j,4], sep="")
    # Possibility of being denovo: 0/0 0/0 1/0
    if(data_pass_site_filter[j,13]==HomRef && data_pass_site_filter[j,19]==HomRef &&
       (data_pass_site_filter[j,7]==HetAlt1 || data_pass_site_filter[j,7]==HetAlt2)){
      denovo[j]= "TRUE"
      type_denovo[j]= "HomRef"
    } else {
      denovo[j]= "FALSE"
      type_denovo[j]= "NA"
    }
    ################################
  } else if (length(alt_nucl)==2){ # We have on ref and 2 alt
    # Set homozygote Alt (2 possibilities)
    HomAlt1 <- paste(alt_nucl[1], "/", alt_nucl[1], sep="")
    HomAlt2 <- paste(alt_nucl[2], "/", alt_nucl[2], sep="")

    # Set heterozygotes
    HetAlt1.1 <- paste(data_pass_site_filter[j,4], "/", alt_nucl[1], sep="")
    HetAlt2.1 <- paste(alt_nucl[1], "/", data_pass_site_filter[j,4], sep="")
    HetAlt1.2 <- paste(data_pass_site_filter[j,4], "/", alt_nucl[2], sep="")
    HetAlt2.2 <- paste(alt_nucl[2], "/", data_pass_site_filter[j,4], sep="")
    HetAlt1.3 <- paste(alt_nucl[1], "/", alt_nucl[2], sep="")
    HetAlt2.3 <- paste(alt_nucl[2], "/", alt_nucl[1], sep="")

    # Possibility of being denovo:
    # 0/0 0/0 1/0
    # 0/0 0/0 2/0
    if(data_pass_site_filter[j,13]==HomRef && data_pass_site_filter[j,19]==HomRef &&
       (data_pass_site_filter[j,7]==HetAlt1.1 || data_pass_site_filter[j,7]==HetAlt2.1 ||
        data_pass_site_filter[j,7]==HetAlt1.2 || data_pass_site_filter[j,7]==HetAlt2.2)){
      denovo[j]= "TRUE"
      type_denovo[j]= "HomRef"
    } else {
      denovo[j]= "FALSE"
      type_denovo[j]= "NA"
    }
    ##################################
  } else { # We have on ref and 3 alt
    # Set homozygote Alt (3 possibilities)
    HomAlt1 <- paste(alt_nucl[1], "/", alt_nucl[1], sep="")
    HomAlt2 <- paste(alt_nucl[2], "/", alt_nucl[2], sep="")
    HomAlt3 <- paste(alt_nucl[3], "/", alt_nucl[3], sep="")

    # Set heterozygotes
    HetAlt1.1 <- paste(data_pass_site_filter[j,4], "/", alt_nucl[1], sep="")
    HetAlt2.1 <- paste(alt_nucl[1], "/", data_pass_site_filter[j,4], sep="")
    HetAlt1.2 <- paste(data_pass_site_filter[j,4], "/", alt_nucl[2], sep="")
    HetAlt2.2 <- paste(alt_nucl[2], "/", data_pass_site_filter[j,4], sep="")
    HetAlt1.3 <- paste(data_pass_site_filter[j,4], "/", alt_nucl[3], sep="")
    HetAlt2.3 <- paste(alt_nucl[3], "/", data_pass_site_filter[j,4], sep="")
    HetAlt1.4 <- paste(alt_nucl[1], "/", alt_nucl[2], sep="")
    HetAlt2.4 <- paste(alt_nucl[2], "/", alt_nucl[1], sep="")
    HetAlt1.5 <- paste(alt_nucl[1], "/", alt_nucl[3], sep="")
    HetAlt2.5 <- paste(alt_nucl[3], "/", alt_nucl[1], sep="")
    HetAlt1.6 <- paste(alt_nucl[2], "/", alt_nucl[3], sep="")
    HetAlt2.6 <- paste(alt_nucl[3], "/", alt_nucl[2], sep="")

    # Possibility of being denovo:
    # 0/0 0/0 1/0
    # 0/0 0/0 2/0
    # 0/0 0/0 3/0
    if(data_pass_site_filter[j,13]==HomRef && data_pass_site_filter[j,19]==HomRef &&
       (data_pass_site_filter[j,7]==HetAlt1.1 || data_pass_site_filter[j,7]==HetAlt2.1 ||
        data_pass_site_filter[j,7]==HetAlt1.2 || data_pass_site_filter[j,7]==HetAlt2.2 ||
        data_pass_site_filter[j,7]==HetAlt1.3 || data_pass_site_filter[j,7]==HetAlt2.3)){
      denovo[j]= "TRUE"
      type_denovo[j]= "HomRef"
    } else {
      denovo[j]= "FALSE"
      type_denovo[j]= "NA"
    }
  }
}


# Format the table
data_pass_site_denovo1 <- data_pass_site_filter[,c(1,2,3)]
data_pass_site_denovo1$denovo <- sapply(denovo, paste0, collapse=",")
data_pass_site_denovo1$typedenovo <- sapply(type_denovo, paste0, collapse=",")
data_pass_site_denovo2 <- data_pass_site_filter[,c(seq(4,length(data_pass_site_filter)))]
data_pass_site_denovo <- cbind(data_pass_site_denovo1, data_pass_site_denovo2)
# Mendelian violation
data_MV <- data_pass_site_denovo[which(data_pass_site_denovo$denovo=="TRUE"),]
########################################################################################

# Apply filters ##########################

# Work only on single SNPs for now:
# Keep only the HomRef and the one SNPs ###############################################
data_pass_HomRef <- data_MV

# Then keep only data that have a good allelic balance  ###############################
AB_filter <- numeric(length = nrow(data_pass_HomRef))
for(i in seq(1, nrow(data_pass_HomRef))){
  idx_other_homo <- 2
  nb_alt <- strsplit(as.character(data_pass_HomRef[i,10]), split=",")[[1]][idx_other_homo]
  tot <- as.numeric(data_pass_HomRef[i,11])
  AB=as.numeric(as.numeric(nb_alt)/tot)
  if(AB=="NaN"){print("Nan")} else{
    if(AB<AB_min || AB>AB_max){
      AB_filter[i] <- "my_ind_filter"
    } else if(AB>=AB_min && AB<=AB_max){
      AB_filter[i] <- "PASS"
    }
  }
}

data_pass_HomRef_AB <- data_pass_HomRef[0,]
for(i in seq(1, nrow(data_pass_HomRef))){
  if(AB_filter[i]=="PASS"){
    data_pass_HomRef_AB <- rbind(data_pass_HomRef_AB, data_pass_HomRef[i,])
  }
}


# Then filter the depth  ############################################################
DP_filter <- numeric(length = nrow(data_pass_HomRef_AB))
for(i in seq(1, nrow(data_pass_HomRef_AB))){
  if(data_pass_HomRef_AB[i,11]<DP_min || data_pass_HomRef_AB[i,11]>DP_max ||
     data_pass_HomRef_AB[i,17]<DP_min || data_pass_HomRef_AB[i,17]>DP_max ||
     data_pass_HomRef_AB[i,23]<DP_min || data_pass_HomRef_AB[i,23]>DP_max){
    DP_filter[i] <- "my_ind_filter"
  } else if(data_pass_HomRef_AB[i,11]>=DP_min && data_pass_HomRef_AB[i,11]<=DP_max &&
            data_pass_HomRef_AB[i,17]>=DP_min && data_pass_HomRef_AB[i,17]<=DP_max &&
            data_pass_HomRef_AB[i,23]>=DP_min && data_pass_HomRef_AB[i,23]<=DP_max){
    DP_filter[i] <- "PASS"}
}


data_pass_HomRef_AB_DP <- data_pass_HomRef_AB[0,]
for(i in seq(1, nrow(data_pass_HomRef_AB))){
  if(DP_filter[i]=="PASS"){
    data_pass_HomRef_AB_DP <- rbind(data_pass_HomRef_AB_DP, data_pass_HomRef_AB[i,])
  }
}


# Finally filter for the genome quality ##############################################
GQ_filter <- numeric(length = nrow(data_pass_HomRef_AB_DP))
for(i in seq(1, nrow(data_pass_HomRef_AB_DP))){
  if(data_pass_HomRef_AB_DP[i,12]<GQ_lim || data_pass_HomRef_AB_DP[i,18]<GQ_lim || data_pass_HomRef_AB_DP[i,24]<GQ_lim){
    GQ_filter[i] <- "my_ind_filter"
  } else if(data_pass_HomRef_AB_DP[i,12]>=GQ_lim && data_pass_HomRef_AB_DP[i,18]>=GQ_lim && data_pass_HomRef_AB_DP[i,24]>=GQ_lim){
    GQ_filter[i] <- "PASS"
    }
}


data_pass_HomRef_AB_DP_GQ <- data_pass_HomRef_AB_DP[0,]
for(i in seq(1, nrow(data_pass_HomRef_AB_DP))){
  if(GQ_filter[i]=="PASS"){
    data_pass_HomRef_AB_DP_GQ <- rbind(data_pass_HomRef_AB_DP_GQ, data_pass_HomRef_AB_DP[i,])
  }
}
##############################################
data_denovo <- data_pass_HomRef_AB_DP_GQ
nb_denovo <- c(name, nrow(data_denovo))
##############################################
write.table(data_denovo, paste0(direct_denovo, "data_denovo_",name,".tab"), sep="\t")
write(nb_denovo,file=paste0(direct_denovo, "denovo.txt"),append=TRUE, ncolumns=2)
## The END ###################################
