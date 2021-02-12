# Script about parental age influence
library(RColorBrewer)

# Import and format data:
data_mut <- read.csv("mutation_rate.txt", sep =" ", header=FALSE) # mutation_rate.txt is obtained during the detection of de novo mutation
age_father=read.csv("male_age.txt", sep =" ", header=FALSE) # one column with father age
age_mother=read.csv("female_age.txt", sep =" ", header=FALSE) # one column with mother age
fa=read.csv("father.txt", sep =" ", header=FALSE) # one column with father ID
mo=read.csv("mother.txt", sep =" ", header=FALSE) # one column with mother age
data=cbind(data_mut[,1:2],age_father, age_mother, fa, mo)

# Parental correlation on mutation rate ######
# Simple correlation:

# Linear model for both age
fit <- lm(data[,2] ~ data[,3] + data[,4], data=data)
summary(fit)

# Predict mu with this model
mu_list=vector()
for(i in seq(1,nrow(data))){
    mu = 1.355e-09+(4.588e-10*data[i,3])+(7.936e-11*data[i,4])
    mu_list=c(mu_list, mu)
}
cor.test(mu_list, data[,2], method="pearson")

fit <- lm(data[,2] ~ mu_list)
summary(fit) # show results

png("model.png", width = 1300, height = 850)
par(mar=c(10,15,4,2), mgp=c(3, 3, 0))
plot(mu_list, data[,2], pch=19, xlab="", ylab="", xlim=c(4e-9,1.2e-8), ylim=c(4e-9,1.2e-8),cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n")
axis(1, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08, 1.2e-08),label=c("0.4","0.5","0.6", "0.7", "0.8","0.9","1.0","1.1","1.2"), cex.axis=4)
axis(2, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08, 1.2e-08),label=c("0.4","0.5","0.6", "0.7", "0.8","0.9","1.0","1.1","1.2"),cex.axis=4, las=2)
mtext(expression("Expected mutation rate x 10"^-8), 1, line=8, cex=3.5)
mtext(expression("Observed mutation rate x 10"^-8), 2, line=9, cex=3.5)
abline(3.765e-13,1.000e+00, lty=2, lwd=3, col="grey40")
text(c(1.15e-08,1.155e-08), c(0.9e-08, 0.95e-08), labels = c("r = 0.54", "p = 0.016"), font=2, cex=3.5)
text(1.00e-8, 1.18e-8, expression(paste(mu, " = 1.36e"^-9, "+ 7.94e"^-11, "× a"[mat], "+ 4.59e"^-10, "× a"[pat])), font=2, cex=2)
dev.off()

# Father correlation ##################################################

linearMod_father <- lm(data[,2] ~ data[,3], data=data)
summary(linearMod_father)
## confidence interval with predict
a <- data[,3]
m <- data[,2]
linearMod_father <- lm(m ~ a)
newx = seq(4,17,by = 0.5)
conf_interval <- predict(linearMod_father, newdata=data.frame(a=newx), interval="confidence",level = 0.95)

png("father_age_pred.png", width = 1300, height = 850)
par(mar=c(9,23,4,2), mgp=c(3, 3, 0))
plot(data[,3], data[,2], pch=19, xlab="", ylab="", xlim=c(3,16), cex=3, cex.lab=2.5, cex.axis=2.5, col="blue", xaxt="n", yaxt="n")
bluetrans <- rgb(204, 229, 255, 255, maxColorValue=255)
polygon(c(newx, rev(newx)), c(conf_interval[,2], rev(conf_interval[,3])),
     col = bluetrans, border = NA)
points(data[,3], data[,2], pch=19, cex=3, cex.lab=2.5, cex.axis=2.5, col="blue")
abline(1.022e-09, 5.393e-10, lwd=3, lty=2, col="darkblue")
axis(1, cex.axis=4)
axis(2, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08),label=c(expression("0.4 x 10"^-8), expression("0.5 x 10"^-8), expression("0.6 x 10"^-8), expression("0.7 x 10"^-8), expression("0.8 x 10"^-8), expression("0.9 x 10"^-8), expression("1.0 x 10"^-8), expression("1.1 x 10"^-8)),cex.axis=4, las=2)
mtext("Paternal age", 1, line=7, cex=4.5)
mtext("Mutation rate (per generation)", 2, line=20, cex=4.5)
text(c(4.1,4.6), c(1.12e-08, 1.06e-08), labels = c(expression(paste(R[adj]^2," = 0.23")), "p = 0.021"), font=2, cex=3.5)
dev.off()

# Mother effect #################################################
linearMod_mother <- lm(data[,2] ~ data[,4], data=data)
summary(linearMod_mother)
a <- data[,4]
m <- data[,2]
linearMod_mother <- lm(m ~ a)
newx = seq(3,17,by = 0.5)
conf_interval <- predict(linearMod_mother, newdata=data.frame(a=newx), interval="confidence",level = 0.95)

png("mother_age_pred.png", width = 1300, height = 850)
par(mar=c(9,23,4,2), mgp=c(3, 3, 0))
plot(data[,4], data[,2], pch=19, xlab="", ylab="", xlim=c(3,16), cex=3, cex.lab=2.5, cex.axis=2.5, col="red2", xaxt="n", yaxt="n")
pinktrans <- rgb(255, 190, 190, 255, maxColorValue=255)
polygon(c(newx, rev(newx)), c(conf_interval[,2], rev(conf_interval[,3])),
     col = pinktrans, border = NA)
points(data[,4], data[,2], pch=19, cex=3, cex.lab=2.5, cex.axis=2.5, col="red2")
abline(6.200e-09, 1.818e-10, lwd=3, lty=2, col="red3")
axis(1, cex.axis=4)
axis(2, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08),label=c(expression("0.4 x 10"^-8), expression("0.5 x 10"^-8), expression("0.6 x 10"^-8), expression("0.7 x 10"^-8), expression("0.8 x 10"^-8), expression("0.9 x 10"^-8), expression("1.0 x 10"^-8), expression("1.1 x 10"^-8)),cex.axis=4, las=2)
mtext("Maternal age", 1, line=7, cex=4.5)
mtext("Mutation rate (per generation)", 2, line=20, cex=4.5)
text(c(4.1,4.6), c(1.12e-08, 1.06e-08), labels = c(expression(paste(R[adj]^2," = 0.09")), "p = 0.111"), font=2, cex=3.5)
dev.off()

## Age parents and phasing ##########################################
poo <- read.csv("assigned_each.tab", sep ="\t") # obtained from the phasing script
poo_p_p <- poo$prop_paternal[-20]
poo_p_m <- poo$prop_maternal[-20]
poo_p <- poo$paternal[-20]
poo_m <- poo$maternal[-20]
poo_n_p <- poo$paternal[-20]/poo$tot_pat[-20]
poo_n_m <- poo$maternal[-20]/poo$tot_mat[-20]
poo_np_p <- (poo_n_p)/(poo_n_p+poo_n_m)
poo_np_m <- (poo_n_m)/(poo_n_p+poo_n_m)
poo_i_p <- round((poo$unassigned+poo$assigned)[-20]*poo_np_p)
poo_i_m <- round((poo$unassigned+poo$assigned)[-20]*poo_np_m)
##
data_p <- cbind(data, poo_p_p, poo_p_m)
colnames(data_p)<-c("name", "mut", "p_age", "m_age", "p", "m", "poo_p_p", "poo_p_m")
lim_inf=min(data_p$p_age, data_p$m_age)
lim_sup=max(data_p$p_age, data_p$m_age)
##

## MODEL on upscaled LINEAR
data_i <- cbind(data, poo_i_p, poo_i_m)
colnames(data_i)<-c("name", "mut", "p_age", "m_age", "p", "m", "poo_i_p", "poo_i_m")
lim_inf=min(data_p$p_age, data_p$m_age)
lim_sup=max(data_p$p_age, data_p$m_age)

reg_age_f <- lm(data_i$poo_i_p ~ data_i$p_age)
summary(reg_age_f)
reg_age_m <- lm(data_i$poo_i_m ~ data_i$m_age)
summary(reg_age_m)
cor.test(data_i$poo_i_p, data_i$p_age, method="pearson")
cor.test(data_i$poo_i_m, data_i$m_age, method="pearson")
confint(reg_age_f)
confint(reg_age_m)

## confidence interval
a1=data_i$p_age
m=data_i$poo_i_p
fit<- lm(m ~ a1)
new1 = seq(4,17,by = 1)
conf_interval <- predict(fit, newdata=data.frame(a1=new1), interval="confidence",level = 0.95)
a_f=new1
conf_f=conf_interval
#
a1=data_i$m_age
m=data_i$poo_i_m
fit<- lm(m ~ a1)
new1 = seq(3,17,by = 1)
conf_interval <- predict(fit, newdata=data.frame(a1=new1), interval="confidence",level = 0.95)
a_m=new1
conf_m=conf_interval

# Plot model
png("infer_given.png", width = 1300, height = 850)
par(mar=c(9,12,4,2), mgp=c(3, 3, 0))
plot(data_i$p_age, data_i$poo_i_p, col="blue", pch=19, ylim=c(0,40), xlim=c(0, lim_sup), xlab="", ylab="",cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n")
bluetrans <- rgb(204, 229, 255, 255, maxColorValue=255)
polygon(c(a_f, rev(a_f)), c(conf_f[,2], rev(conf_f[,3])), col = bluetrans, border = NA)
pinktrans <- rgb(255, 190, 190, 255, maxColorValue=255)
polygon(c(a_m, rev(a_m)), c(conf_m[,2], rev(conf_m[,3])), col = pinktrans, border = NA)
points(data_i$p_age, data_i$poo_i_p, col="blue", pch=19,cex=3)
points(data_i$m_age, data_i$poo_i_m, col="red2", pch=19,cex=3)
segments(4, 4.8399 + (1.8364*4), 20, 4.8399 + (1.8364*20), lwd=3, lty=2, col="blue")
text(c(14.5,15), c(21.5, 18.5), labels = c(expression(paste(R[adj]^2," = 0.41")), "p = 0.002"), font=2, cex=3, col="blue")
segments(3, 4.6497+ (0.3042*3), 20, 4.6497+ (0.3042*20), lwd=3, lty=2, col="red2")
text(c(14.5,15), c(4.5, 1.5), labels = c(expression(paste(R[adj]^2," = -0.01")), "p = 0.38"), font=2, cex=3, col="red2")
arrows(3,20,3,4.6497+ (0.3042*3), lwd=2, lty=1, col="red2")
arrows(4,20,4,4.8399 + (1.8364*4), lwd=2, lty=1, col="blue")
text(3.5, 22, "Puberty", cex=3.5, font=2)
legend("topleft", legend=c("Maternal", "Paternal"), pch=c(19,19), col=c("red2", "blue"), text.font=2, cex=3.5, bty = "n")
axis(1, cex.axis=4)
axis(2, cex.axis=4, las=2)
mtext("Parental age", 1, line=7, cex=4)
mtext("Upscaled phased mutations", 2, line=8, cex=4)
dev.off()


## Test of the model
# Predict mu with this model
nb_list=vector()
for(i in seq(1,nrow(data))){
    nb = (4.8399 + 1.8364*data[i,3]) + (4.6497+ 0.3042*data[i,4])
    nb_list=c(nb_list, nb)
}
call <- read.csv("../../de_novo_mutation/callability.txt", sep =" ", header=FALSE)
callability=mean(call[,2])
fnr=0.04020232294960846
mu_infer=(nb_list*(1-0.10887096774193548))/(2*callability*(1-fnr))

cor.test(mu_infer, data[,2], method="pearson")
fit <- lm(data[,2] ~ mu_infer)
summary(fit) # show results

png("model_upscaled_mu.png", width = 1300, height = 850)
par(mar=c(10,15,4,2), mgp=c(3, 3, 0))
plot(mu_infer, data[,2], pch=19, xlab="", ylab="", xlim=c(4.0e-9,1.2e-8), ylim=c(4.0e-9,1.2e-8),cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n")
axis(1, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08, 1.2e-08),label=c("0.4", "0.5", "0.6", "0.7", "0.8","0.9","1.0","1.1","1.2"), cex.axis=4)
axis(2, at=c(4.0e-09, 5.0e-09, 6.0e-09, 7.0e-09, 8.0e-09, 9.0e-09, 1.0e-08, 1.1e-08, 1.2e-08),label=c("0.4", "0.5", "0.6", "0.7", "0.8","0.9","1.0","1.1","1.2"),cex.axis=4, las=2)
mtext(expression("Expected mutation rate x 10"^-8), 1, line=8, cex=3.5)
mtext(expression("Observed mutation rate x 10"^-8), 2, line=9, cex=3.5)
abline(-1.064e-09, 1.277, lty=2, lwd=3, col="grey40")
text(c(1.15e-08,1.155e-08), c(1.0e-08, 1.05e-08), labels = c("r = 0.54", "p = 0.016"), font=2, cex=3.5)
text(0.77e-8, 1.165e-8, expression(paste(mu, " =")), font=2, cex=2)
segments(0.8e-8, 1.165e-8, 1.2e-8, 1.165e-8, lwd=2, lty=1)
text(1.0e-8, 1.19e-8, expression(paste("(4.65 + 0.30 × a"[mat], " + 4.84 + 1.84 × a"[pat],") x (1-0.089)")), font=2, cex=2)
text(1.0e-8, 1.14e-8, expression(paste("2 x 2351302179 x (1-0.0402)")), font=2, cex=2)
dev.off()


### Same with Poisson #################################
## MODEL on upscaled POISSON
data_i <- cbind(data, poo_i_p, poo_i_m)
colnames(data_i)<-c("name", "mut", "p_age", "m_age", "p", "m", "poo_i_p", "poo_i_m")
lim_inf=min(data_p$p_age, data_p$m_age)
lim_sup=max(data_p$p_age, data_p$m_age)

# Regression
reg_age_f <- glm(data_i$poo_i_p ~ data_i$p_age, family = poisson)
summary(reg_age_f)
##reg_age_f <- glm(data_i$poo_i_p ~ data_i$p_age*data_i$p, family = poisson)
reg_age_m <- glm(data_i$poo_i_m ~ data_i$m_age, family = poisson)
summary(reg_age_m)
##
a1=data_i$p_age
m=data_i$poo_i_p
fit<- glm(m ~ a1, family = poisson)
new1 = seq(4,17,by = 0.5)
conf_interval <- predict(fit, newdata=data.frame(a1=new1), interval="confidence", type="response", level = 0.95, se.fit = TRUE)
a_f=new1
conf_f=conf_interval$fit
uci_f=conf_interval$fit+1.96*conf_interval$se.fit
lci_f=conf_interval$fit-1.96*conf_interval$se.fit
#
# Female
a1=data_i$m_age
m=data_i$poo_i_m
fit<- glm(m ~ a1, family = poisson)
new1 = seq(3,17,by = 0.5)
conf_interval <- predict(fit, newdata=data.frame(a1=new1), interval="confidence", type="response", level = 0.95, se.fit = TRUE)
a_m=new1
conf_m=conf_interval$fit
uci_m=conf_interval$fit+1.96*conf_interval$se.fit
lci_m=conf_interval$fit-1.96*conf_interval$se.fit
#

##plot
png("infer_given_poisson.png", width = 1300, height = 850)
par(mar=c(9,12,4,2), mgp=c(3, 3, 0))
plot(data_i$p_age, data_i$poo_i_p, col="blue", pch=19, ylim=c(0,40), xlim=c(0, lim_sup), xlab="", ylab="",cex=3, cex.lab=2.5, cex.axis=2.5, xaxt="n", yaxt="n")
bluetrans <- rgb(204, 229, 255, 255, maxColorValue=255)
polygon(c(a_f, rev(a_f)), c(lci_f, rev(uci_f)),col = bluetrans, border = NA)
pinktrans <- rgb(255, 190, 190, 255, maxColorValue=255)
polygon(c(a_m, rev(a_m)), c(lci_m, rev(uci_m)),col = pinktrans, border = NA)
lines(a_f, conf_f, lwd=3, lty=2, col="blue")
lines(a_m, conf_m, lwd=3, lty=2, col="red2")
points(data_i$p_age, data_i$poo_i_p, col="blue", pch=19, cex=3)
points(data_i$m_age, data_i$poo_i_m, col="red2", pch=19, cex=3)
##text(c(9,6.2), c(39, 36), labels = c("nb = exp(2.14+0.089xage)", "p < 0.001"), font=2, cex=3, col="blue")
##text(c(13,10.8), c(6, 3), labels = c("nb = exp(1.35+0.067xage)", "p < 0.01"), font=2, cex=3, col="red2")
legend("topleft", legend=c("Maternal", "Paternal"), pch=c(19,19), col=c("red2", "blue"), text.font=2, cex=3.5, bty = "n")
axis(1, cex.axis=4)
axis(2, cex.axis=4, las=2)
mtext("Parental age", 1, line=7, cex=4)
mtext("Upscaled phased mutations", 2, line=8, cex=4)
dev.off()


### Box plot contribution ##################################################################
png("norm_prop_box.png", width = 1300, height = 850)
##png("norm_prop_box.png", width = 950, height = 850)
par(mar=c(9,17,4,2), mgp=c(3, 3, 0))
boxplot(poo_np_m, poo_np_p, ylim=c(0,1), names=c("",""), ylab="", cex=3, cex.lab=2.5, cex.axis=2.5, col=c("red2","blue"), yaxt="n")
axis(2, cex.axis=3.5, las=2)
mtext("Proportion of phased", 2, line=12, cex=4)
mtext(expression(italic(de~novo)~mutations), 2, line=9, cex=4)
axis(1, at=c(1,2), label=c("Maternal", "Paternal"), cex.axis=4)
dev.off()
#

# Stats
se=sd(poo_np_p)/sqrt(length(poo_np_p))
se
mean(poo_np_p)-(2.1*se)
mean(poo_np_p)+(2.1*se)

# t test ####################################################
t.test(poo_np_p,poo_np_m)
