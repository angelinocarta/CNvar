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
diver<-read.csv("data/diversificationALL.csv", header=T) #this is a diversification file created by myself but only used to extract the clades information
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
diver<-read.csv("data/medusafeatures.csv", header=T)
diver2<-unique(merge(diver, ho, by="genus")[,c(1,8,22)])
diver2<-aggregate(r.mean ~ family, diver2, mean)

data<-merge(data, diver2, by="family")

diver3<-read.csv("data/family.stem.ageSB.csv", header=T)
data<-merge(data, diver3, by="family")
data$age<-data$stem.age

data$div.rate0<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0)
data$div.rate0.9<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0.9)
data<-subset(data, div.rate0>0) ##### removing families with a single species
data<-subset(data, r.mean>0.01) ##### removing outliers: families with very low diversification
plot(log(data$SD2), log(data$r.mean))

#####
ttt<-read.tree("data/Family_SBmcc.v2.tre")##this is the SANTIAGO tree
cont<-read.csv("data/order.constraint.csv", header=T)
data <- data[match(ttt$tip.label, data$family), ]
data.cont<-merge(data, cont, by="order")
data.cont <- data.cont[order(data.cont$ordine.APG),]
t<-drop.tip(ttt, as.character(data.cont$family))
ttt<-drop.tip(ttt, as.character(t$tip.label))
ttt$node.label<-paste("N",1:3241,sep="")
ttt<-force.ultrametric(ttt, method="nnls")
ttt<-ladderize(ttt, right = TRUE)
ttt<- read.tree(text = write.tree((ttt)))
ttt<-ape::rotate(ttt, 335)
ttt<- read.tree(text = write.tree((ttt)))
ttt

#####
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
tree<-ttt
p <- ggtree(tree, layout="fan", open.angle=10, size=0.5, ladderize = F, right = T  )#, branch.length = 'none')

grp <- list(asterids     =  subset(data, clade=="asterids")$family,
            monocots     =  subset(data, clade=="monocots")$family,
            otherAngiosperms     =  subset(data, clade=="other Angiosperms")$family,
            rosids     =  subset(data, clade=="rosids")$family)

p2<-groupOTU(p, grp) + aes(color=group) +
  scale_color_manual(values=c("orange", "olivedrab3", "blue", "purple")) +
  theme(legend.position="right")

# Then add a bar layer outside of the tree.
p3 <- p2 + 
  new_scale_fill() +
  geom_fruit(
    data=data,
    geom=geom_bar,
    mapping=aes(y=family, x=SD, fill=clade),  #The 'Abundance' of 'dat1' will be mapped to x
    pwidth=0.4,
    stat="identity",
    orientation="y", # the orientation of axis.
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      hjust=0,  # adjust the horizontal position of text of axis.
      text.size = 2,
      title = "Number of cytotypes"
    ),
    grid.params=list() # add the grid line of the external bar plot.
  ) + 
  scale_fill_manual(
    values=c("orange", "olivedrab3", "blue", "purple"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=6)
  ) +
  theme(#legend.position=c(0.96, 0.5), # the position of legend.
    legend.background=element_rect(fill=NA), # the background of legend.
    legend.title=element_text(size=7), # the title size of legend.
    legend.text=element_text(size=6), # the text size of legend.
    legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
  )

pdf("results/Fig.1_circle.pdf", width=10, height=10)
p3
dev.off()
