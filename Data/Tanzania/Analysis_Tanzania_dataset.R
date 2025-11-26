library(viridis)

setwd('~/Desktop/MiNI_PW/Xenstats/Data')

file<-'data_tanzania.csv'

data.tanzania<-read.csv(file)

colnames(data.tanzania)

n<-dim(data.tanzania)[1]

############################################
#### PLOTTING QUANTITATIVE PHENOTYPES ######
############################################

par(mfrow=c(1,2))

boxplot(data.tanzania$Height.cm~data.tanzania$Ethnic.Group,xlab='Ethnic group',ylab='Height (in cm)',las=1,col=viridis(4))

boxplot(data.tanzania$Temperature~data.tanzania$Ethnic.Group,xlab='Ethnic group',ylab='Body temperature (in Celsius)',las=1,col=viridis(4))

###########################################################
#### TESTING ASSOCIATION BETWEEN rs1801033 AND HEIGHT ######
###########################################################

genotypes<-unique(data.tanzania$rs1801033)

table(data.tanzania$rs1801033,useNA='always')

par(mfrow=c(1,1))

boxplot(data.tanzania$Height.cm~data.tanzania$rs1801033,xlab='Genotype',ylab='Height (in cm)',las=1,col=viridis(4),main='rs1801033')

data.tanzania$additive.effects.rs1801033<-rep(NA,n)

data.tanzania$additive.effects.rs1801033[data.tanzania$rs1801033=='AA']<-0

data.tanzania$additive.effects.rs1801033[data.tanzania$rs1801033=='AC']<-1

data.tanzania$additive.effects.rs1801033[data.tanzania$rs1801033=='CC']<-2

table(data.tanzania$additive.effects.rs180133,useNA='always')

fit.m1<-lm(Height.cm~additive.effects.rs180133,data=data.tanzania)

#### Score test ###

res<-summary(fit.m1)

coef(res)

p.value.score.test<-coef(res)['additive.effects.rs180133','Pr(>|t|)']

#### Likelihood ratio rest ###

loglik.m1<-logLik(fit.m1)

fit.m0<-lm(Height.cm~1,data=data.tanzania[is.na(data.tanzania$additive.effects.rs180133)==F,])

loglik.m0<-logLik(fit.m0)

stat<-(-2)*(loglik.m0-loglik.m1)

p.value.lr.test<-pchisq(stat,df=1,lower.tail=F)

###########################################################
#### TESTING ASSOCIATION BETWEEN rs1801033 AND HEIGHT ######
###########################################################

fit.m1.extended<-lm(Height.cm~additive.effects.rs180133+Ethnic.Group,data=data.tanzania)

summary(fit.m1.extended)

loglik.m1.extended<-logLik(fit.m1.extended)

fit.m0.extended<-lm(Height.cm~ Ethnic.Group,data=data.tanzania[is.na(data.tanzania$additive.effects.rs180133)==F,])

loglik.m0.extended<-logLik(fit.m0.extended)

stat<-(-2)*(loglik.m0.extended-loglik.m1.extended)

p.value.lr.test<-pchisq(stat,df=1,lower.tail=F)

###########################################################
#### TESTING ASSOCIATION BETWEEN rs6874639 AND ANAEMIA ####
###########################################################

genotypes<-unique(data.tanzania$rs6874639)

table(data.tanzania$rs6874639,useNA='always')

par(mfrow=c(1,1))

tabela1<-table(data.tanzania$rs6874639,data.tanzania$Anaemia)

proportion.per.genotype<-tabela1[,2]/rowSums(tabela1)

barplot(proportion.per.genotype,beside=T,ylab='Height (in cm)',xlab='genotypes',las=1,col=viridis(4),main='rs6874639')

data.tanzania$additive.effects.rs6874639<-rep(NA,n)

data.tanzania$additive.effects.rs6874639[data.tanzania$rs6874639=='AA']<-0

data.tanzania$additive.effects.rs6874639[data.tanzania$rs6874639=='AG']<-1

data.tanzania$additive.effects.rs6874639[data.tanzania$rs6874639=='GG']<-2

data.tanzania$Anaemia.bin<-ifelse(data.tanzania$Anaemia=='Yes',1,ifelse(is.na(data.tanzania$Anaemia)==F,0,NA))

fit.m1.probit<-glm(Anaemia.bin~ additive.effects.rs6874639,data=data.tanzania,family=binomial(link='probit'))

fit.m1.logit<-glm(Anaemia.bin~ additive.effects.rs6874639,data=data.tanzania,family=binomial(link='logit'))

res.probit<-summary(fit.m1.probit)

p.value.score.test.probit<-coef(res.probit)['additive.effects.rs6874639','Pr(>|z|)']

res.logit<-summary(fit.m1.logit)

p.value.score.test.logit<-coef(res.logit)['additive.effects.rs6874639','Pr(>|z|)']

coef(fit.m1.probit)

coef(fit.m1.logit)/1.70

###########################################################
#### TESTING ASSOCIATION BETWEEN rs3024500 AND ANAEMIA ####
###########################################################

genotypes<-unique(data.tanzania$rs3024500)

table(data.tanzania$rs3024500,useNA='always')

par(mfrow=c(1,1))

tabela1<-table(data.tanzania$rs3024500,data.tanzania$Anaemia)

proportion.per.genotype<-tabela1[,2]/rowSums(tabela1)

barplot(proportion.per.genotype,beside=T,ylab='Height (in cm)',xlab='genotypes',las=1,col=viridis(4),main='rs3024500')

data.tanzania$additive.effects.rs3024500<-rep(NA,n)

data.tanzania$additive.effects.rs3024500[data.tanzania$rs3024500=='AA']<-0

data.tanzania$additive.effects.rs3024500[data.tanzania$rs3024500=='AG']<-1

data.tanzania$additive.effects.rs3024500[data.tanzania$rs3024500=='GG']<-2

fit.m1.probit<-glm(Anaemia.bin~ additive.effects.rs3024500,data=data.tanzania,family=binomial(link='probit'))

fit.m1.logit<-glm(Anaemia.bin~ additive.effects.rs3024500,data=data.tanzania,family=binomial(link='logit'))

res.probit<-summary(fit.m1.probit)

p.value.score.test.probit<-coef(res.probit)['additive.effects.rs3024500','Pr(>|z|)']

res.logit<-summary(fit.m1.logit)

p.value.score.test.logit<-coef(res.logit)['additive.effects.rs3024500','Pr(>|z|)']

coef(fit.m1.probit)

coef(fit.m1.logit)/1.70

###########################################################
#### TESTING HWE: rs1801033, rs6874639, and rs3024500  ####
###########################################################

allele.freq<-sum(data.tanzania$additive.effects.rs1801033,na.rm=T)/(2*sum(is.na(data.tanzania$additive.effects.rs1801033)==F))

res.rs1801033<-chisq.test(table(data.tanzania$additive.effects.rs1801033),p=c((1-allele.freq)^2,2*allele.freq*(1-allele.freq),allele.freq^2))

res.rs1801033

rbind(res.rs1801033$observed,res.rs1801033$expected)

allele.freq<-sum(data.tanzania$additive.effects.rs6874639,na.rm=T)/(2*sum(is.na(data.tanzania$additive.effects.rs6874639)==F))

res.rs6874639<-chisq.test(table(data.tanzania$additive.effects.rs6874639),p=c((1-allele.freq)^2,2*allele.freq*(1-allele.freq),allele.freq^2))

res.rs6874639

rbind(res.rs6874639$observed,res.rs6874639$expected)

allele.freq<-sum(data.tanzania$additive.effects.rs3024500,na.rm=T)/(2*sum(is.na(data.tanzania$additive.effects.rs3024500)==F))

res.rs3024500 <-chisq.test(table(data.tanzania$additive.effects.rs3024500),p=c((1-allele.freq)^2,2*allele.freq*(1-allele.freq),allele.freq^2))

res.rs3024500

rbind(res.rs3024500 $observed,res.rs3024500$expected)
