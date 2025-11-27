###################################################################
# Packages needed                                             #####
# MASS needed for correcting p-values                         #####
# betareg needed for beta regression                          #####
# If you don't have the packages installed                    #####
# in your computer,                                           #####
# Use install.packages(c('MASS','betareg'), dependencies=T)   #####
###################################################################

library(MASS)

library(betareg)

####################################################
########### LOADING DATA.         ##################
####################################################

setwd('~/Desktop/MiNI_PW/Xenstats/Data')

file.data<-'data_ace_ace2.csv'

data<-read.csv(file.data,header=T)

boxplot(data[,'cg23232263']~data[,'Disease'],ylim=c(0,1),las=1)

table(data$Study,data$Disease)

#### Selecting data from Herrera et al ####

data<-data[data$Study==4,]

colnames(data)

data[1,]

dim(data)

### 109 individuals and 31 columns 

######## Loading annotation ########

file.annotation<-'ace_ace2_annotation.csv'

data.annotation<-read.csv(file.annotation,header=T)

colnames(data.annotation)

data.annotation[1,]

## TSS - transcription start site

####################################################
########### Analysis of cg18877734 (beta values) ###
########### Mann-Whitney / T test                ###
####################################################

par(bg='white')

boxplot(data[,'cg18877734']~data[,'Disease'],ylim=c(0.9,1),las=1,xlab='Disease status',ylab='Proportion of methylation (beta values)',main='cg18877734',col=c('light green','salmon'))

wilcox.test(data[,'cg18877734']~data[,'Disease'])$p.value

t.test(data[,'cg18877734']~data[,'Disease'])$p.value	

res1<-lm(cg18877734~Disease+Female,data=data)

summary(res1)

####################################################
########### Analysis of cg18877734 (M values)    ###
########### Mann-Whitney / T test                ###
####################################################

data$M.values.cg18877734<-log(data[,'cg18877734']/(1-data[,'cg18877734']))

par(bg='white')

boxplot(data$M.values.cg18877734~data$Disease,ylim=c(2,5),las=1,xlab='Disease status',ylab='log odds (M values)',main='cg18877734',col=c('light green','salmon'))

wilcox.test(data[,'M.values.cg18877734']~data[,'Disease'])$p.value

t.test(data[,'M.values.cg18877734']~data[,'Disease'])$p.value	

res2<-lm(M.values.cg18877734~Disease+Female,data=data)

summary(res2)


####################################################
########### Analysis of cg18877734 (beta values) ###
########### Beta regression                      ###
####################################################

res3<-betareg(cg18877734~Disease+Female,data=data)	

summary(res3)