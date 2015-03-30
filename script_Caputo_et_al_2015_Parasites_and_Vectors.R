# packages ####
library(plyr)
library(ggplot2)
library(geoR)
library(SoDA)
library(gstat)
library(glmmADMB)

# MET          ####

# Import Data  ####
Dati2 <- read.csv("Emerging_R.csv", sep=";",header=TRUE,na.strings="na")
summary(Dati2)
# removing NA
dati_nona <- na.omit(Dati2)
summary(dati_nona)

# date
levels(dati_nona$data_raccolta)
# data preparation
dati_nona$settimana    <- factor(dati_nona$settimana)
levels(dati_nona$sito) <- c("UN","S","W","E","N")
#IGIENE - W
#POLITICHE -N
#FARMACIA - S
#ORTO - E

# exploratory analysis
boxplot(albo_fem~dentro_fuori,dati_nona)
boxplot(pip_fem~dentro_fuori,dati_nona)
boxplot(albo_fem~dentro_fuori,dati_nona[dati_nona$sito=="UN",])
boxplot(pip_fem~dentro_fuori,dati_nona[dati_nona$sito=="UN",])
# total inside/outside catches by sex and species
sum(dati_nona$albo_fem)
sum(dati_nona$albo_mas)
sum(dati_nona$pip_fem)
sum(dati_nona$pip_mas)


# GLMM outside ####
OUT <- droplevels(subset(dati_nona, dentro_fuori=="O"))
summary(OUT)
table(OUT$sito)
table(OUT$settimana)
tapply(OUT$albo_fem,OUT$sito,mean)


oun <- droplevels(subset(OUT,sito=="UN"))
sum(oun$albo_fem,oun$albo_mas)/sum(oun$pip_fem,oun$pip_mas)

boxplot(albo_fem~sito, dati_nona)

# mean outside catches by site

media <- ddply(OUT,c("sito"),summarize,
               albofem=round(mean(albo_fem),1),seae = round(sd(albo_fem)/sqrt(length(albo_fem)),1),
               pipfem=round(mean(pip_fem),1),  ster = round(sd(pip_fem)/sqrt(length(pip_fem)),1),
               albmale=round(mean(albo_mas),1), seam = round(sd(albo_mas)/sqrt(length(albo_mas)),1),
               pipmale=round(mean(pip_mas),1),  sepm = round(sd(pip_mas)/sqrt(length(pip_mas)),1)
               )
media
media[c(5,4,2,3,1),]

# total outside catches by sex and species
all<-ddply(OUT,c("dentro_fuori"),summarize,alb=sum(albo_fem),pip=sum(pip_fem),
           albmale=sum(albo_mas),pipmale=sum(pip_mas))
all

OUT$tot_albo <- OUT$albo_fem + OUT$albo_mas
OUT$tot_pip  <- OUT$pip_fem  + OUT$pip_mas

change_reference<-function(y, x){
  library(glmmADMB)
  OUT$sito <- relevel(OUT$sito,x)
  a <- substitute(y~sito)
  mm.out <- glmmadmb(a,random=~1|settimana, family="nbinom",data=OUT)
  par(mfrow=c(2,2))
  plot(residuals(mm.out,type="pearson"),ylab="Pearson Residuals");abline(0,0)
  plot(fitted(mm.out),residuals(mm.out),ylab="Pearson Residuals",xlab="Fitted values")
  plot(OUT$sito,residuals(mm.out,type="pearson"),ylab="Pearson Residuals",xlab="Sito")
  plot(OUT$settimana,residuals(mm.out),ylab="Pearson Residuals",xlab="Settimana")
  return(summary(mm.out))}


# ________Albopictus  Outside ######
change_reference(albo_fem,"UN")
change_reference(albo_fem,"S") # 
change_reference(albo_fem,"W") #  
change_reference(albo_fem,"E") # 
change_reference(albo_fem,"N") #  

change_reference(albo_mas,"UN")
change_reference(albo_mas,"S") 
change_reference(albo_mas,"W") 
change_reference(albo_mas,"E")
change_reference(albo_mas,"N")  

# ________Pipiens  Outside ######
change_reference(pip_fem,"UN")
change_reference(pip_fem,"S") # 
change_reference(pip_fem,"W") #   
change_reference(pip_fem,"E") #  
change_reference(pip_fem,"N") #  

change_reference(pip_mas,"UN")
change_reference(pip_mas,"S") #  
change_reference(pip_mas,"W") #   
change_reference(pip_mas,"E") #  
change_reference(pip_mas,"N") # 
# n politiche
# e orto 
# w igiene
# s farmicia



# Fisher test inside   #####

In <- droplevels(subset(dati_nona, dentro_fuori=="I"))
summary(In)
iun<-droplevels(subset(In,sito=="UN"))
sum(iun$albo_fem,iun$albo_mas)/sum(iun$pip_fem,iun$pip_mas)

boxplot(albo_fem~sito,In)

table(In$settimana)
table(In$sito)

# mean inside catches by site
media <- ddply(In,c("sito"),summarize,
               albofem=round(mean(albo_fem),1),seae = round(sd(albo_fem)/sqrt(length(albo_fem)),1),
               pipfem=round(mean(pip_fem),1),  ster = round(sd(pip_fem)/sqrt(length(pip_fem)),1),
               albmale=round(mean(albo_mas),1), seam = round(sd(albo_mas)/sqrt(length(albo_mas)),1),
               pipmale=round(mean(pip_mas),1),  sepm = round(sd(pip_mas)/sqrt(length(pip_mas)),1),
               ae=round(mean(albo_fem+albo_mas),1),  se = round(sd(albo_fem+albo_mas)/sqrt(length(albo_fem+albo_mas)),1),
               pip=round(mean(pip_fem+pip_mas),1),  sepip = round(sd(pip_fem+pip_mas)/sqrt(length(pip_fem+pip_mas)),1)
               )
media

# build presence/absence data from inner side met
In$presenza_fem  <- ifelse(In$albo_fem > 0,1,0)
In$presenza_mas  <- ifelse(In$albo_mas > 0,1,0)
In$presenza_femc <- ifelse(In$pip_fem > 0,1,0)
In$presenza_masc <- ifelse(In$pip_mas > 0,1,0)
prespip          <- In$presenza_masc + In$presenza_femc
presalbo         <- In$presenza_mas  + In$presenza_fem
In$prespip       <- ifelse(prespip  > 0,1,0)
In$presalbo      <- ifelse(presalbo > 0,1,0)

table(droplevels(subset(In,sito %in% c("UN","N") ))$sito)

table(In$sito,In$presenza_fem)
table(In$sito,In$presenza_mas)
table(In$sito,In$presenza_femc)
table(In$sito,In$presenza_masc)

In$Riclas_Location <- factor(ifelse(In$sito=="UN","UN","T"))
table(In$Riclas_Location,In$presenza_fem)
table(In$Riclas_Location,In$presenza_mas)
table(In$Riclas_Location,In$presalbo)
table(In$Riclas_Location,In$prespip)


In_AO <- droplevels(subset(In, sito=="UN" | sito=="E"))
fisher.test(table(In_AO$sito, In_AO$presenza_fem))
fisher.test(table(In_AO$sito, In_AO$presenza_mas))
fisher.test(table(In_AO$sito, In_AO$presenza_femc))
fisher.test(table(In_AO$sito, In_AO$presenza_masc))


In_AA<-droplevels(subset(In, Riclas_Location =="UN" | Riclas_Location =="T"))
fisher.test(table(In_AA$presenza_fem,In_AA$Riclas_Location))
fisher.test(table(In_AA$presenza_mas,In_AA$Riclas_Location))
fisher.test(table(In_AA$presenza_femc,In_AA$Riclas_Location))
fisher.test(table(In_AA$presenza_masc,In_AA$Riclas_Location))


In_OP<-droplevels(subset(In, sito=="UN" | sito=="N"))
fisher.test(table(In_OP$sito, In_OP$presenza_fem))
fisher.test(table(In_OP$sito, In_OP$presenza_mas))
fisher.test(table(In_OP$sito, In_OP$presenza_femc))
fisher.test(table(In_OP$sito, In_OP$presenza_masc))


In_AF<-droplevels(subset(In, sito=="UN" | sito=="W"))
fisher.test(table(In_AF$sito, In_AF$presenza_fem))
fisher.test(table(In_AF$sito, In_AF$presenza_mas))
fisher.test(table(In_AF$sito, In_AF$presenza_femc))
fisher.test(table(In_AF$sito, In_AF$presenza_masc))


In_AF<-droplevels(subset(In, sito=="UN" | sito=="S"))
fisher.test(table(In_AF$sito, In_AF$presenza_fem))
fisher.test(table(In_AF$sito, In_AF$presenza_mas))
fisher.test(table(In_AF$sito, In_AF$presenza_femc))
fisher.test(table(In_AF$sito, In_AF$presenza_masc))


# Adulticide ####

# Import Data  ####
Dati <- read.csv("SAP_AN_R.csv", sep=";",header=TRUE,na.strings="na")

# tolgo NA
Dati_NoNA <- Dati[-c(6, 64, 268, 478),]
str(Dati_NoNA)
Dati_NoNA$Settimana <- factor(Dati_NoNA$Settimana)
Dati_NoNA$griglia[is.na(Dati_NoNA$griglia)] <- 0
Dati_NoNA$griglia <- factor(Dati_NoNA$griglia)


xy_all           <- geoXY(Dati_NoNA$lat, Dati_NoNA$long, unit=1)
Dati_NoNA[30:31] <- xy_all
names(Dati_NoNA)[30:31] <- c("Lat","Long")

# GLM NB relationship between ST- and CBT #####

dbase0 <-ddply(Dati_NoNA,c("griglia","Settimana"),transform)
head(dbase0)
dbase  <- droplevels(subset(dbase0, Settimana %in% c(3,4,6,7,9,10)))
dbase2 <- data.frame(yST  = dbase$albo_fem[dbase$trappola=="STs"],
    			 xCBT = dbase$albo_fem[dbase$trappola=="Tombino"],
					 sito = dbase$sito[dbase$trappola=="Tombino"])


levels(dbase2$sito) <- c("Untreated","Treated")
dbase2$sito<-relevel(dbase2$sito,"Treated")
mod1 <- glm.nb(yST ~ xCBT*sito, data = dbase2)
summary(mod1)


modtreat <- glm.nb(yST ~ xCBT, data = droplevels(subset(dbase2, sito=="Treated")))
summary(modtreat)
anova(modtreat, test="Chisq")

moduntreat <- glm.nb(yST ~ xCBT, data = droplevels(subset(dbase2, sito=="Untreated")))
summary(moduntreat)
anova(moduntreat, test="Chisq")


dbase  <- droplevels(subset(dbase0, Settimana %in% c(3,4,6,7,9,10)))
names(dbase)

dbase3 <- data.frame(yST  = dbase$pip_fem[dbase$trappola=="STs"],
					 xCBT = dbase$pip_fem[dbase$trappola=="Tombino"],
					 sito = dbase$sito[dbase$trappola=="Tombino"])

modpip1 <- glm.nb(yST ~ xCBT*sito, data = dbase3)
summary(modpip1)
moduntpip <- glm.nb(yST ~ xCBT, data = droplevels(subset(dbase3, sito=="anatomia")))
summary(moduntpip)

modtreatpip <- glm.nb(yST ~ xCBT, data = droplevels(subset(dbase3, sito=="sapienza")))
summary(modtreatpip)
anova(modtreatpip,test="Chisq")

# GLMMs NB, effect of treatments and collection method ######

#mosquito counts 
ddply(Dati_NoNA,c("sito","trappola"),summarize,
      af = round(mean(albo_fem),1),s1 =round(sd(albo_fem)/sqrt(length(albo_fem)),1),
      cf = round(mean(pip_fem),1), s3 =round(sd(pip_fem)/sqrt(length(pip_fem)),1),
      am = round(mean(albo_mas),1),s2 =round(sd(albo_mas)/sqrt(length(albo_mas)),1),
      cm = round(mean(pip_mas),1), s4 =round(sd(pip_mas)/sqrt(length(pip_mas)),1)
      )

# glmm ####
albo_femglmm<-glmmadmb(albo_fem~sito*trappola+(1|Settimana), 
                 family="nbinom",
                 data= Dati_NoNA)
summary(albo_femglmm)


(sum(Dati_NoNA$albo_fem[Dati_NoNA$sito=="anatomia"])-
  sum(Dati_NoNA$albo_fem[Dati_NoNA$sito=="sapienza"]))/sum(Dati_NoNA$albo_fem[Dati_NoNA$sito=="anatomia"])
(exp(2.3379)-exp( 2.3379-0.9173))/exp(2.3379)
(exp(2.3379)-exp( 2.3379-0.8641))/exp(2.3379)

(exp(2.3379)-exp( 2.3379-0.9173))/exp(2.3379) # difference betwen un an t in st
(exp(2.3379-0.8641)-exp( 2.3379-0.9173-0.8641))/exp(2.3379-0.8641) # difference betwen un an t in cbt

dati_glmm <- data.frame(R=resid(albo_femglmm,type="pearson"), X= Dati_NoNA$Lat, Y= Dati_NoNA$Long,set=Dati_NoNA$Settimana,tr=Dati_NoNA$trappola,sito=Dati_NoNA$sito)
coordinates(dati_glmm)<-c("X","Y")
variogramma.t<-variogram(R ~ 1,dati_glmm)
plot(variogramma.t,col="black",pch=20,cex=2)

plot(fitted(albo_femglmm),resid(albo_femglmm,type="pearson"))
plot(resid(albo_femglmm,type="pearson"))

pip_femglmm<-glmmadmb(pip_fem~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati_NoNA)
summary(pip_femglmm)


(exp(0.822)-exp( 0.822-0.353))/exp(0.822)
(exp( 0.822+0.475)-exp(0.822))/exp(0.822+0.475)
(exp( 0.822+0.475)-exp(0.822))/exp(0.822)

(exp( 0.822)-exp(0.822-0.353 ))/exp(0.822) # st fem unt vs treat
(exp( 0.822+0.475)-exp(0.822-0.353+0.475-0.802 ))/exp(0.822+0.475) # cbt fem unt vs treat



albo_masglmm<-glmmadmb(albo_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati_NoNA)
summary(albo_masglmm)
(exp(1.7224)-exp( 1.7224-0.9692))/exp(1.7224)
(exp(1.7224)-exp( 1.7224-0.6552))/exp(1.7224)

(exp(1.7224-0.6552)-exp( 1.7224-0.9692-0.6552))/exp(1.7224-0.6552)

pip_masglmm<-glmmadmb(pip_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati_NoNA)
summary(pip_masglmm)
(exp(0.0602)-exp( 0.0602-0.2418))/exp(0.0602)
(exp( 0.0602+0.5898)-exp(0.0602))/exp(0.0602+0.5898)

(exp( 0.0602)-exp(0.0602-0.2418))/exp(0.0602)
(exp( 0.0602+0.5898)-exp(0.0602+0.5898-0.8418))/exp(0.0602+0.5898)

# Henderson ####     
media <- ddply(Dati_NoNA,c("sito","Settimana"),summarize, 
               MeanAe = round(mean(albo_fem+albo_mas,na.rm=T),1),
               seAe   = round(sd(albo_fem+albo_mas,na.rm=T)/sqrt(length(albo_fem+albo_mas)),1),
               MeanCx = round(mean(pip_fem+pip_mas,na.rm=T),1),
               seCx   = round(sd(pip_fem+pip_mas,na.rm=T)/sqrt(length(pip_fem+pip_mas)),1))
media
round(media[,3:5],2)

length(Dati_NoNA[Dati_NoNA$Settimana==2 & Dati_NoNA$sito=="sapienza",]$albo_fem)

tr = media[13,3]/media[12,3]
un = media[3,3]/media[2,3]

percentcontrol1 = 100 - ((tr/un)*100)
round(percentcontrol1,1)

tr2 = media[17,3]/media[16,3]
un2 = media[7,3]/media[6,3]

percentcontrol2 = 100 - ((tr2/un2)*100)
round(percentcontrol2,1)


tr = media[13,5]/media[12,5]
un = media[3,5]/media[2,5]

percentcontrol1 = 100 - ((tr/un)*100)
round(percentcontrol1,1)

tr2 = media[17,5]/media[16,5]
un2 = media[7,5]/media[6,5]

percentcontrol2 = 100 - ((tr2/un2)*100)
round(percentcontrol2,1)


cbtred <- droplevels(subset(Dati_NoNA,trappola=="Tombino"))
stred <- droplevels(subset(Dati_NoNA,trappola=="STs"))
media <- ddply(cbtred,c("sito","Settimana"),summarize, 
               MeanAe = round(mean(albo_fem+albo_mas,na.rm=T),1),
               seAe   = round(sd(albo_fem+albo_mas,na.rm=T)/sqrt(length(albo_fem+albo_mas)),1),
               MeanCx = round(mean(pip_fem+pip_mas,na.rm=T),1),
               seCx   = round(sd(pip_fem+pip_mas,na.rm=T)/sqrt(length(pip_fem+pip_mas)),1))
media
round(media[,3:8],2)
round(media[c(2,3,6,7),c(3,5)],2) # untreat
round(media[c(12,13,16,17),c(3,5)],2) # treat

length(Dati_NoNA[Dati_NoNA$Settimana==2 & Dati_NoNA$sito=="sapienza",]$albo_fem)

tr = media[13,3]/media[12,3]
un = media[3,3]/media[2,3]

percentcontrol1 = 100 - ((tr/un)*100)
round(percentcontrol1,1)

tr2 = media[17,3]/media[16,3]
un2 = media[7,3]/media[6,3]

percentcontrol2 = 100 - ((tr2/un2)*100)
round(percentcontrol2,1)


tr = media[13,5]/media[12,5]
un = media[3,5]/media[2,5]

percentcontrol1 = 100 - ((tr/un)*100)
round(percentcontrol1,1)

tr2 = media[17,5]/media[16,5]
un2 = media[7,5]/media[6,5]

percentcontrol2 = 100 - ((tr2/un2)*100)
round(percentcontrol2,1)

posttr <- subset(Dati_NoNA, (Settimana == 3 | Settimana ==7) & sito == "sapienza")
mean(posttr$albo_fem,na.rm=T)

pretr <- subset(Dati_NoNA, (Settimana == 2 | Settimana ==6) & sito == "sapienza")
mean(pretr$albo_fem,na.rm=T)

trtot <- mean(posttr$albo_fem,na.rm=T)/mean(pretr$albo_fem,na.rm=T)


postun <- subset(Dati_NoNA, (Settimana == 3 | Settimana ==7) & sito == "anatomia")
mean(postun$albo_fem,na.rm=T)

preun <- subset(Dati_NoNA, (Settimana == 2 | Settimana ==6) & sito == "anatomia")
mean(preun$albo_fem,na.rm=T)

untot <- mean(postun$albo_fem,na.rm=T)/mean(preun$albo_fem,na.rm=T)

percentcontrol = 100 - ((trtot/untot)*100)
percentcontrol

(27.2+81.7)/2
#where T is the post application mean divided by the pre application mean in the treatment site 
#U is the post application mean divided by the pre application mean in the control 


# GLMMs Grid #####
sap<-droplevels(subset(Dati_NoNA,sito=="sapienza"))
albo_femgriglia<-glmmadmb(albo_fem~trappola*griglia,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= sap)
summary(albo_femgriglia)

dati_glmm<-data.frame(R=resid(albo_femgriglia,type="pearson"), X= sap$Lat, Y= sap$Long,set=sap$Settimana,tr=sap$trappola)
coordinates(dati_glmm)<-c("X","Y")
variogramma.t<-variogram(R ~ 1,dati_glmm)
plot(variogramma.t,col="black",pch=20,cex=2)
plot(resid(albo_femgriglia,type="pearson"))
plot(fitted(albo_femgriglia),resid(albo_femgriglia,type="pearson"))
plot(sap$Long,resid(albo_femgriglia,type="pearson"))
plot(sap$Lat,resid(albo_femgriglia,type="pearson"))
plot(sap$griglia,resid(albo_femgriglia,type="pearson"))

pip_femgriglia<-glmmadmb(pip_fem~trappola*griglia,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= sap)
summary(pip_femgriglia)


albo_masgriglia<-glmmadmb(albo_mas~trappola*griglia,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= sap)
summary(albo_masgriglia)


pip_masgriglia<-glmmadmb(pip_mas~trappola*griglia,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= sap)
summary(pip_masgriglia)

# temporal effect of adulticide treatmentsE #####
Dati<-droplevels(subset(Dati_NoNA, as.numeric(Settimana) <7))
albofem7<-glmmadmb(albo_fem~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati)
summary(albofem7)

pipfem7<-glmmadmb(pip_fem~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati)
summary(pipfem7)

albomas7<-glmmadmb(albo_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati)
summary(albomas7)

pipmas7<-glmmadmb(pip_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati)
summary(pipmas7)


Dati3<-droplevels(subset(Dati_NoNA, as.numeric(Settimana) <3))
albofem3<-glmmadmb(albo_fem~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati3)
summary(albofem3)

pipfem3<-glmmadmb(pip_fem~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati3)
summary(pipfem3)

albomas3<-glmmadmb(albo_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati3)
summary(albomas3)

pipmas3<-glmmadmb(pip_mas~sito*trappola,
                 random=~1|Settimana, 
                 family="nbinom",
                 zeroInflation=F,
                 data= Dati3)
summary(pipmas3)


## ratio ######
num<-ddply(Dati_NoNA,c("sito","trappola"),summarize,taf=sum(albo_fem),tam=sum(albo_mas),tpf=sum(pip_fem),tpm=sum(pip_mas))
num
PerformanceAN <-
  matrix(c(num$taf[1],num$taf[2],num$tpf[1],num$tpf[2]),
         nrow = 2,
         dimnames = list("Trap" = c("ST", "CBT"),
                         "Species" = c("Alb", "Pip")))
PerformanceSAP <-
  matrix(c(num$taf[3],num$taf[4],num$tpf[3],num$tpf[4]),
         nrow = 2,
         dimnames = list("Trap" = c("ST", "CBT"),
                         "Species" = c("Alb", "Pip")))
chisq.test(PerformanceAN)
chisq.test(PerformanceSAP)
fisher.test(PerformanceAN)
fisher.test(PerformanceSAP)
prop.test(PerformanceAN)
prop.test(PerformanceSAP)

Dati_NoNA$uno<-rep(1,nrow(Dati_NoNA))
numsd<-ddply(Dati_NoNA,c("sito","trappola"),summarize,
             taf=mean(pip_fem),
             tafd=(sd(pip_fem)/sqrt(sum(uno))))
round(numsd[,3:4],2)


num2<-ddply(Dati_NoNA,c("sito","trappola"),summarize,
            taf=mean(albo_fem),tam=mean(albo_mas),tpf=mean(pip_fem),tpm=mean(pip_mas))
num2
round(num2[,3:6],2)
round(num2$taf/num2$tpf,2)
round(num2$tam/num2$tpm,2)
num2

femmCBTan<-c(num2$taf[2],num2$tpf[2])
binom.test(femmCBTan,0.5)
chisq.test(femmCBTan)

femmCBTan<-c(num2$tam[2],num2$tpm[2])
binom.test(femmCBTan,0.5)
chisq.test(femmCBTan)
chisq.test(cbind(c(1094,432),c(231,369)))
chisq.test(cbind(c(1094,432),c(231,369)))



