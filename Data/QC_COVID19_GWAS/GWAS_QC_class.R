library(MASS)

setwd('~/Desktop/COVID19/Before_filtering/')

##############################################################################
####################### CHECK ALLELE FREQUENCIES #############################
##############################################################################

file<-'plink.frq'

data.freq<-read.table(file,header=T,stringsAsFactors=F)

colnames(data.freq)

dim(data.freq)

summary(100*data.freq$MAF)

sum(data.freq$MAF==0)

sum(data.freq$MAF<0.01)

sum(data.freq$MAF<0.025)

hist(100*data.freq$MAF,xlim=c(0,50),las=1,xlab='Minor allele frequency (%)',ylab='proportion',main='',lwd=2,axes=T,freq=F,nclass=40,ylim=c(0,0.15))

plot(ecdf(100*data.freq$MAF),xlim=c(0,50),las=1,xlab='Minor allele frequency (%)',ylab='empirical cumulative distribution',main='',lwd=2,axes=T)

#########################################################################
####################### Missing data ####################################
#######################     SNP      ####################################
#########################################################################

file<-'plink.lmiss'

data.missing.snp<-read.table(file,header=T,stringsAsFactors=F)

summary(data.missing.snp$'F_MISS')

ecdf.miss<-ecdf(100*data.missing.snp$'F_MISS')

plot(knots(ecdf.miss),log10(1-ecdf.miss(knots(ecdf.miss))),type='l',xlim=c(0,20),las=1,xlab='Missing data frequency (%)',ylab='1-ecdf',main='SNP',lwd=2,axes=T,ylim=c(-6,0))

#########################################################################
####################### Missing data ####################################
#######################  Individual  ####################################
#########################################################################

file<-'plink.imiss'

data.missing.ind<-read.table(file,header=T,stringsAsFactors=F)

summary(data.missing.ind$'F_MISS')

ecdf.miss.ind<-ecdf(100*data.missing.ind$'F_MISS')

plot(knots(ecdf.miss.ind),log10(1-ecdf.miss.ind(knots(ecdf.miss.ind))),type='s',xlim=c(0,7),las=1,xlab='Missing data frequency (%)',ylab='1-ecdf',main='Individual',lwd=2,axes=T,ylim=c(-4,0))

#########################################################################
########################## CHECK HWE ####################################
#########################################################################

data.hwe.all<-read.table('plink.hwe2.all',header=T,stringsAsFactors=F)

colnames(data.hwe)

plot(ecdf(data.hwe.all$P),las=1,xlab='p-value (HWE test)',ylab='Cumulative probability distribution',main='All individuals',lwd=2)

#abline(c(0,1),col='red',lwd=3)

lines(c(0,1),c(0,1),col='red',lwd=3)

p.values.bonf<-p.adjust(data.hwe.all$P,method='bonferroni')

sum(p.values.bonf<0.05)

sum(P.values.fdr<0.05)

#########################################################################
####################### CHECK HETEROZYGOSITY ############################
#########################################################################

data.heterozygosity<-data.hwe.all[,'O.HET.']

summary(data.heterozygosity)

boxplot(data.heterozygosity,las=1,ylab='heterozygosity',ylim=c(0,0.6),col='mistyrose1')

#########################################################################
########################## CHECK SEX ####################################
#########################################################################

file<-'plink.sexcheck'

data.sexcheck<-read.table(file,header=T,stringsAsFactors=F)

table(data.sexcheck$PEDSEX,data.sexcheck$SNPSEX)
