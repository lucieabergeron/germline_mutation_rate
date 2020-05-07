# Import data ##########################################################################
data <- read.csv(paste0(direct_denovo,"output.table.fnr.", name), sep ="\t")
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
# Keep only the filtered ones
data_pass_site_filter <- data_order[which(data_order$FILTER=="PASS"),]
########################################################################################

# Check for allelic balances  ##########################################################
AB_filter <- numeric(length = nrow(data_pass_site_filter))
for(i in seq(1, nrow(data_pass_site_filter))){
  idx_other_homo <- 2
  nb_tot <- data_pass_site_filter[i,9]
  nb_alt_min <- min(strsplit(as.character(data_pass_site_filter[i,8]), split=",")[[1]])
  AB=as.integer(nb_alt_min)/nb_tot
  if(AB=="NaN"){print("Nan")} else{
    if(AB<AB_min || AB>AB_max){
      AB_filter[i] <- "my_ind_filter"
    } else if(AB>=AB_min && AB<=AB_max){
      AB_filter[i] <- "PASS"
    }
  }
}

########################################################################################
den=nrow(data_pass_site_filter)
nom=length(which(AB_filter=="my_ind_filter"))
fnr=nom/den
stat=c(name, den, nom, fnr)

write(stat,file=paste0(direct_denovo, "fnr.txt"),append=TRUE, ncolumns=4)

## The END ###################################
