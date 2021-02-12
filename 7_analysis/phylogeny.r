# Script for phylogeny dating

# 1. Ne estimation ################################################################################################
data <- read.csv("mutation_rate.txt", sep =" ", header=FALSE)

average=mean(data[,2])
average
ne= 0.00247 /(4*mu)

# 2. Dating phylogeny ################################################################################################
library(RColorBrewer)
library(dplyr)
library(ape)
library(phytools)

gt_m = 11
# Distances
d_Pongo_abelii <- 0.0162551430672998
d_Homo_sapiens <- 0.00577676246248343
d_Pan_troglodytes <- 0.00612226407385174
d_h_c <- 0.00928793113723326
d_ape <- 0.0111854267754341
d_Macaca_mulatta <-0.00260875810745442
d_Macaca_f <- 0.0022282726931979
d_macaca_macaca <-0.00464525362885664
d_baboon <- 0.00686281100087727
d_maca_bab <- 0.00399650084545229
d_Chlorocebus_sabaeus <- 0.0111626595169656
d_m_g <- 0.0246717581129659

# Using only the new estimated on the branches in question
mu_rhesus= ((4.8399 + 1.8364 *12+ 4.6497+ 0.3042 *10)*(1-0.10887096774193548)/(2*2351302179*(1-0.04020232294960846)))
mu_m = (mu_rhesus/11)*1000000
mu_rhesus_min= ((25.158457+5.065827)*(1-0.10887096774193548)/(2*2351302179*(1-0.04020232294960846)))
CI_min= (mu_rhesus_min/11)*1000000
mu_rhesus_max= ((28.59583+10.317069)*(1-0.10887096774193548)/(2*2351302179*(1-0.04020232294960846)))
CI_max= (mu_rhesus_max/11)*1000000

# Macaca-macaca
t_m = (d_Macaca_mulatta)/(mu_m)
t_m_min = (d_Macaca_mulatta)/(CI_min)
t_m_max = (d_Macaca_mulatta)/(CI_max)
# Macaca-baboon
t_m_b = (d_Macaca_mulatta+d_macaca_macaca)/(mu_m)
t_m_b_min = (d_Macaca_mulatta+d_macaca_macaca)/(CI_min)
t_m_b_max = (d_Macaca_mulatta+d_macaca_macaca)/(CI_max)
# Macaca-green
t_m_g = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab)/(mu_m)
t_m_g_min = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab)/(CI_min)
t_m_g_max = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab)/(CI_max)

# Macaca with great apes
t_m_a = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab+d_m_g)/(mu_m)
t_m_a_min = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab+d_m_g)/(CI_min)
t_m_a_max = (d_Macaca_mulatta+d_macaca_macaca+d_maca_bab+d_m_g)/(CI_max)

##
# Babon rate:
mu_b = 0.549e-9*1000000
d_baboon/mu_b
# Green rate:
mu_g = 1.11e-9*1000000
d_Chlorocebus_sabaeus/mu_g
# Homo:
mu_h= 0.43e-9*1000000
(d_ape+d_h_c+d_Homo_sapiens)/mu_h
# Chimp
mu_c= 0.64e-9*1000000
(d_ape+d_h_c+d_Pan_troglodytes)/mu_c

# Model
mu_nwm=2.7e-9*1000000
(d_Macaca_mulatta+d_macaca_macaca)/mu_m
(d_maca_bab+d_m_g)/mu_g
d_maca_bab/mu_g
d_m_g/mu_nwm

####
tree_t_new_reg<-read.tree("tree_timed_rhesus_estimate_last_reg.nwk")
tree_t_topo<-read.tree("tree_timed_rhesus_estimate_topo_reg.nwk")
tree_t_sp<-read.tree("tree_timed_rhesus_estimate_sp_reg.nwk")

CI= rbind(c(t_m_a_min,t_m_a_max),
          c(t_m_g_min,t_m_g_max),
          c(t_m_b_min,t_m_b_max),
          c(t_m_min,t_m_max))

# speciation
minus_m=(2*79874*11)/1000000
minus_cat=(2*355000*11)/1000000

sp_t_m = t_m -minus_m
sp_t_m_a =t_m_a-minus_cat

#
png("phylo_diff_sp_n_reg_last.png", width = 1200, height = 750)
par(xaxt="n",yaxt="n",mar=c(5.1,6.1,2.1,2), mgp=c(5, 1.5, 0))
plotTree.errorbars(tree_t_new_reg, CI, fsize=3, ftype="off", lwd=10, bar.width=0, cex=2.5, xlab="",ylab="", at=seq(0,45,by=5))
#
polygon(c(sp_t_m,sp_t_m,t_m,t_m), c(1,2,2,1), col = "gray80", border="gray80")
polygon(c(sp_t_m_a,sp_t_m_a,t_m_a,t_m_a), c(3.13,5,5,3.13), col = "gray80", border="gray80")
#
axis(1, at=seq(0,55,by=5), cex.axis=2)
#
abline(v=60, lty=2, lwd=4, col="gray65")
abline(v=55, lty=2, lwd=4, col="gray65")
abline(v=50, lty=2, lwd=4, col="gray65")
abline(v=45, lty=2, lwd=4, col="gray65")
abline(v=40, lty=2, lwd=4, col="gray65")
abline(v=35, lty=2, lwd=4, col="gray65")
abline(v=30, lty=2, lwd=4, col="gray65")
abline(v=25, lty=2, lwd=4, col="gray65")
abline(v=20, lty=2, lwd=4, col="gray65")
abline(v=15, lty=2, lwd=4, col="gray65")
abline(v=10, lty=2, lwd=4, col="gray65")
abline(v=5, lty=2, lwd=4, col="gray65")
abline(v=0, lty=2, lwd=4, col="gray65")
#
lines(c(50.09,50.09),c(5,3.15), lwd=10, col="gray52")
lines(c(2.45,2.45),c(1,2), lwd=10, col="gray52")
points(c(50.09),c(4.06), lwd=15,pch=20, col="navy")
points(c(2.45),c(1.5), lwd=15,pch=20, col="navy")
par(new=TRUE)
par(xaxt="n",yaxt="n",mar=c(5.1,5.1,2.1,2), mgp=c(5, 1.5, 0))
## ftype= "reg" "i" "b"
plotTree.errorbars(tree_t_new_reg, CI, fsize=3, ftype="i", lwd=10, bar.width=20, cex=2.5, xlab="",ylab="", at=seq(0,45,by=5),bar.col="darkgreen")
##plotTree.errorbars(tree_t_sp, CI, fsize=3, ftype="i", lwd=10, bar.width=20, cex=2.5, xlab="",ylab="", at=seq(0,45,by=5),bar.col="darkblue")
##CI_s=rbind(cbind(50.09,50.09),c(18.13,18.13),c(11.69,11.69),c(2.45,2.45))
##plotTree.errorbars(tree_t_sp, CI_s, fsize=3, ftype="i", lwd=10,bar.width=20, cex=2.5, xlab="",ylab="", at=seq(0,45,by=5),bar.col="darkblue")
#
par(xaxt="s",yaxt="s",font.lab=4)
axis(1, at=seq(0,60,by=5), cex.axis=2)
#
text(c(10,19.5,30,58.2), c(1.35,2.1,2.98,3.915),labels = c("Macaca", "Papionini", "Cercopithecidae", "Catarrhini"), font=1, cex=2.3)
text(c(10,58.2), c(1.2,3.765),labels = c(expression("(N"[e]*" = 79,874)"), expression("(N"[e]*" = 355,000)")), font=2, cex=1.6)
text(c(-0.1,-2), c(1.62,1.48),labels = c(expression("T"[d]*"=4.20"),expression("T"[s]*"=2.45")), col=c("darkgreen","navy"), font=2, cex=1.8)
text(c(6.5), c(2.4),labels = c(expression("T"[d]*"=11.69")), col=c("darkgreen"), font=2, cex=1.8)
text(c(13), c(3.28),labels = c(expression("T"[d]*"=18.13")), col=c("darkgreen"), font=2, cex=1.8)
text(c(52.5,45.5), c(4.17,3.9),labels = c(expression("T"[d]*"=57.90"),expression("T"[s]*"=50.09")), col=c("darkgreen","navy"), font=2, cex=1.8)
#
# add fossils:
points(c(22,22,33.5), c(3.12,5,4.05), pch=4, lwd=5, col="red")
text(c(28,22,33.5), c(3.22,5.1,4.15),labels = c("Victoriapithecus", "Proconsul", "Aegyptopithecus"), , col=c("red"), font=3, cex=1.8)
dev.off()

# 3. compare estimate with literature ######################################################################################################################
div_me <- c(4.20,11.69,18.13,57.90)
div_lit_nucl <- c(3.53, 8.13, 11.5, 31.6)
div_lit_mito <- c(3.44, 12.17, 14.09, 32.12)

png("estimation_dif_met.png", width = 1300, height = 850)
par(mar=c(9,15,4,2), mgp=c(3, 3, 0))
plot(div_lit_mito, div_me, type='o', pch=19, lwd=4, xlim=c(0,60), ylim=c(0,60), xlab="", ylab="",
     cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n", lty=1)
lines(div_lit_nucl, div_me, type='o', pch=19, lwd=4, xlab="", ylab="", cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n", lty=1, col="gray60")
abline(0,1, lwd=2, lty=1, col="black")
points(div_lit_mito, div_me, pch=21, lwd=3, cex=4, cex.lab=2.5, cex.axis=2.5, col=1, bg=c("aquamarine4", "firebrick1", "dodgerblue4", "darkgoldenrod2"))
points(div_lit_nucl, div_me, pch=21, lwd=3, cex=4, cex.lab=2.5, cex.axis=2.5, col="gray60", bg=c("aquamarine4", "firebrick1", "dodgerblue4", "darkgoldenrod2"))
text(c(7.5, 18, 6.5, 38),c(1, 11.5, 26.5, 56), labels=c("Macaca", "Papionini", "Cercopithecidae","Catarrhini"),
     font=c(4,2,2,2), cex=3, col=c("aquamarine4", "firebrick1", "dodgerblue4", "darkgoldenrod2"))
axis(1, cex.axis=4)
axis(2, cex.axis=4, las=2)
mtext("Time (in Mya) from the molecular clock", 1, line=7, cex=4)
mtext("Time (in Mya) from our estimate", 2, line=9, cex=4)
legend("topleft", c(expression("T"[divergence]*" compare to mitochondrial data"), expression("T"[divergence]*" compare to nuclear data")), col=c("black", "gray60"), lty=c(1,1), lwd=3, cex=2, bty = "n")
dev.off()


