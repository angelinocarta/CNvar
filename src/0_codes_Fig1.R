####### FAMILY LEVEL
# Set the working directory
library(chromer)
library("CCDBcurator")
library(taxonlookup)
library(phytools)
library(geiger)
library(caper)
library(nlme)

data(angiorecordsclean2018)
ho<-add_higher_order()
ho2<-ho[,c(3,4)]
angiorecordsclean2018$genus<- angiorecordsclean2018$Genus
merged1<-merge(angiorecordsclean2018, ho2, by="genus")
merged1$species<-paste(merged1$Genus, merged1$Species, sep="_")

#unique chromosome by family
fam1<- unique(merged1[,c(9,8)])
species.cn<- na.omit(fam1[order(fam1$family, fam1$csome.number),])
species.cn<-aggregate(csome.number ~ family, species.cn, length)
family.chromo.sd.2018<-species.cn
names(family.chromo.sd.2018)<-c("family", "SD")

##add order and clade
diver<-read.csv("diversificationALL.csv", header=T) #this is a diversification file created by myself but only used to extract the clades information
clades<-unique(diver[,c(6,7,8)])
merged<-merge(family.chromo.sd.2018, clades, by="family")

##add SR from TPL and sampling species
sampled<-unique(merged1[,c(9,10)])
sampled<-data.frame(table(sampled$family)) #count the number of sampled species per family
names(sampled)<-c("family", "sampled")
merged<-merge(merged, sampled, by="family")
#merged$SR.davies<-merged$SR
merged$SR<-NULL
richness<-aggregate(number.of.accepted.species~family, ho, sum) #count the number of species per family
names(richness)<-c("family", "SR")
merged<-merge(merged, richness, by="family")

data<-merged
data$sampling.effort<-data$sampled/data$SR
data$sampling.effort[data$sampling.effort > 1] <- 1

data$SD2<-data$SD # diversity of cytotypes

##add MEDUSA per family
diver<-read.csv("medusafeatures.csv", header=T)
diver2<-unique(merge(diver, ho, by="genus")[,c(1,8,22)])
diver2<-aggregate(r.mean ~ family, diver2, mean)

data<-merge(data, diver2, by="family")

diver3<-read.csv("family.stem.ageSB.csv", header=T)
data<-merge(data, diver3, by="family")
data$age<-data$stem.age

data$div.rate0<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0)
data$div.rate0.9<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0.9)
data<-subset(data, div.rate0>0) ##### removing families with a single species
#data$SD2<- data$SD2 +0.0001
#data$div.rate0<- data$div.rate0 +0.0001
#data$div.rate0.9<- data$div.rate0 +0.0001
#data$r.mean<- data$r.mean +0.0001
data<-subset(data, r.mean>0.01) ##### removing outliers: families with very low diversification
plot(log(data$SD2), log(data$r.mean))

#####
ttt<-read.tree("Family_SBmcc.v2.tre")##this is the SANTIAGO tree
cont<-read.csv("order.constraint.csv", header=T)
data <- data[match(ttt$tip.label, data$family), ]
data.cont<-merge(data, cont, by="order")
data.cont <- data.cont[order(data.cont$ordine.APG),]
t<-drop.tip(ttt, as.character(data.cont$family))
ttt<-drop.tip(ttt, as.character(t$tip.label))
ttt$node.label<-paste("N",1:3241,sep="")
ttt<-force.ultrametric(ttt, method="nnls")
ttt<-ladderize(ttt, right = TRUE)
ttt<- read.tree(text = write.tree((ttt)))
ttt<-rotate(ttt, 335)
ttt<- read.tree(text = write.tree((ttt)))
