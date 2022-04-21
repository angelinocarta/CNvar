####### GENUS LEVEL
library(chromer)
library("CCDBcurator")
library(taxonlookup)
library(phytools)
library(geiger)
library(caper)
library(nlme)
library(rr2)

data(angiorecordsclean2018)
angiorecordsclean2018$genus<- angiorecordsclean2018$Genus
angiorecordsclean2018$species<-paste(angiorecordsclean2018$Genus, angiorecordsclean2018$Species, sep="_")

#unique chromosome by genus
genus1<- unique(angiorecordsclean2018[,c(1,7)])
species.cn<- na.omit(genus1[order(genus1$Genus, genus1$csome.number),])
species.cn<-aggregate(csome.number ~ Genus, species.cn, length)
genus.chromo.sd.2018<-species.cn
names(genus.chromo.sd.2018)<-c("genus", "SD")

clades<-read.csv("data/clades.csv", header=T)
merged<-merge(genus.chromo.sd.2018, clades, by="genus")

ho<-add_higher_order()
SR<-ho[,c(3,1)]
names(SR)<-c("genus", "SR")
merged<-merge(merged, SR, by="genus")

sampled<-unique(angiorecordsclean2018[,c(8,9)])
sampled<-data.frame(table(sampled$genus))
names(sampled)<-c("genus", "sampled")
merged<-merge(merged, sampled, by="genus")

data<-merged
data$sampling.effort<-data$sampled/data$SR
data$sampling.effort[data$sampling.effort > 1] <- 1

data$SD2<-data$SD

diver<-read.csv("data/medusafeatures.csv", header=T)
diver2<-aggregate(r.mean ~ genus, diver, mean)
merged<-merge(data, diver2, by="genus")

diver3<-read.csv("data/genera.stem.ageSB.csv", header=T)
merged<-merge(merged, diver3, by="genus")

tree<-read.tree("data/RC_complete_MCCv_2_GENERA.tre") ##SANTIAGO tree
t<-drop.tip(tree, as.character(merged$genus))
ttt<-drop.tip(tree, as.character(t$tip.label))
ttt$node.label<-paste("N",1:ttt$Nnode,sep="")

data <- merged[match(ttt$tip.label, merged$genus), ]
data$age<-data$stem.age
data$div.rate0<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0)
data$div.rate0.9<- bd.ms(time= data$age, n=data$SR, crown=F, epsilon=0.9)
data<-subset(data,div.rate0>0)
data$age<-data$stem.age
plot(log10(data$SD2), log10(data$div.rate0))

######
data$species2<-data$genus
da<-subset(data, clade=="monocots")
dp<-subset(data, clade=="rosids")
dw<-subset(data, clade=="asterids")

cdat <- comparative.data(data=data, phy=ttt, names.col="species2")
signal.test<-phylosig(cdat[[1]], setNames(log(cdat[[2]]$SD2), rownames(cdat[[2]])), method="lambda",test=TRUE)
##
ca <- comparative.data(data=da, phy=ttt, names.col="species2")
cp <- comparative.data(data=dp, phy=ttt, names.col="species2")
cw <- comparative.data(data=dw, phy=ttt, names.col="species2")

r1 <- gls(log10(div.rate0) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cdat[[1]], fixed = FALSE), data = cdat[[2]], method = "ML")
r2 <- gls(log10(div.rate0.9) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cdat[[1]], fixed = FALSE), data = cdat[[2]], method = "ML")
r3 <- gls(log10(r.mean) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cdat[[1]], fixed = FALSE), data = cdat[[2]], method = "ML")
r4 <- gls(log10(SR) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cdat[[1]], fixed = FALSE), data = cdat[[2]], method = "ML")
r5 <- gls(log10(age) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cdat[[1]], fixed = FALSE), data = cdat[[2]], method = "ML")

r1a <- gls(log10(div.rate0) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = ca[[1]], fixed = FALSE), data = ca[[2]], method = "ML")
r2a <- gls(log10(div.rate0.9) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = ca[[1]], fixed = FALSE), data = ca[[2]], method = "ML")
r3a <- gls(log10(r.mean) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = ca[[1]], fixed = FALSE), data = ca[[2]], method = "ML")
r4a <- gls(log10(SR) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = ca[[1]], fixed = FALSE), data = ca[[2]], method = "ML")
r5a <- gls(log10(age) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = ca[[1]], fixed = FALSE), data = ca[[2]], method = "ML")

r1p <- gls(log10(div.rate0) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cp[[1]], fixed = FALSE), data = cp[[2]], method = "ML")
r2p <- gls(log10(div.rate0.9) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cp[[1]], fixed = FALSE), data = cp[[2]], method = "ML")
r3p <- gls(log10(r.mean) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cp[[1]], fixed = FALSE), data = cp[[2]], method = "ML")
r4p <- gls(log10(SR) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cp[[1]], fixed = FALSE), data = cp[[2]], method = "ML")
r5p <- gls(log10(age) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cp[[1]], fixed = FALSE), data = cp[[2]], method = "ML")

r1w <- gls(log10(div.rate0) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cw[[1]], fixed = FALSE), data = cw[[2]], method = "ML")
r2w <- gls(log10(div.rate0.9) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cw[[1]], fixed = FALSE), data = cw[[2]], method = "ML")
r3w <- gls(log10(r.mean) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cw[[1]], fixed = FALSE), data = cw[[2]], method = "ML")
r4w <- gls(log10(SR) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cw[[1]], fixed = FALSE), data = cw[[2]], method = "ML")
r5w <- gls(log10(age) ~ log10(SD2), weights = ~1/sampling.effort, correlation = corPagel(1, phy = cw[[1]], fixed = FALSE), data = cw[[2]], method = "ML")

pdf("results/2.pgls.genera.SB.pdf", width=4, height=8)
par(mfrow=c(5,1))
par(mar = c(0,4,0,4), oma = c(5, 2, 3, 2))
plot(jitter(cdat[[2]]$SD2, factor=0,amount=0.05), (cdat[[2]]$div.rate0), pch = 21, bg="gray", col="darkgray", cex = 1.1 * sqrt(cdat[[2]]$sampling.effort), xlab = "", ylab = "", xlim=c(1,250), ylim=c(0.001, 1), yaxt="n", xaxt="n", log="xy")
legend(45, 0.01, legend=c("all angiosperms", "monocots", "rosids", "asterids"),col=c("blue", "olivedrab3", "purple", "orange"), lty=c(1,2,2,2), cex=0.6, box.lty=0)
at.y <- outer(1:9, 10^(-3:1))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1 , cex.axis=0.75)
abline(r1, col = "blue", lwd=2)
abline(r1a, col = "olivedrab3", lwd=1 , lty=2)
abline(r1p, col = "purple", lwd=1, lty=2)
abline(r1w, col = "orange", lwd=1, lty=2)
plot( jitter(cdat[[2]]$SD2, factor=0,amount=0.05), (cdat[[2]]$div.rate0.9), pch = 21, bg="gray", col="darkgray", cex = 1.1 * sqrt(cdat[[2]]$sampling.effort), xlab = "", ylab = "", xlim=c(1,250), ylim=c(0.0001, 1), yaxt="n", xaxt="n", log="xy")
at.y <- outer(1:9, 10^(-4:1))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1, cex.axis=0.75)
abline(r2, col = "blue", lwd=2)
abline(r2a, col = "olivedrab3", lwd=1, lty=2)
abline(r2p, col = "purple", lwd=1, lty=2)
abline(r2w, col = "orange", lwd=1, lty=2)
plot( jitter(cdat[[2]]$SD2, factor=0,amount=0.05), (cdat[[2]]$r.mean), pch = 21, bg="gray", col="darkgray", cex = 1.1 * sqrt(cdat[[2]]$sampling.effort), xlab = "", ylab = "", xlim=c(1,250), ylim=c(0.01, 2), yaxt="n", xaxt="n", log="xy")
at.y <- outer(1:9, 10^(-2:1))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1, cex.axis=0.75)
abline(r3, col = "blue", lwd=2)
abline(r3a, col = "olivedrab3", lwd=1, lty=2)
abline(r3p, col = "purple", lwd=1, lty=2)
abline(r3w, col = "orange", lwd=1, lty=2)
plot( jitter(cdat[[2]]$SD2, factor=0,amount=0.05), (cdat[[2]]$SR), pch = 21, bg="gray", col="darkgray", cex = 1.1 * sqrt(cdat[[2]]$sampling.effort), xlab = "", ylab = "", xlim=c(1,250), ylim=c(1, 2000), yaxt="n", xaxt="n", log="xy")
at.y <- outer(1:9, 10^(0:5))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1, cex.axis=0.75)
abline(r4, col = "blue", lwd=2)
abline(r4a, col = "olivedrab3", lwd=1, lty=2)
abline(r4p, col = "purple", lwd=1, lty=2)
abline(r4w, col = "orange", lwd=1, lty=2)
plot( jitter(cdat[[2]]$SD2, factor=0,amount=0.05), (cdat[[2]]$age), pch = 21, bg="gray", col="darkgray", cex = 1.1 * sqrt(cdat[[2]]$sampling.effort), xlab = "", ylab = "", xlim=c(1,250), ylim=c(10, 200), yaxt="n", xaxt="n", log="xy")
at.x <- outer(1:9, 10^(0:3))
lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
axis(1, at=at.x, labels=lab.x, las=1, cex.axis=0.75)
at.y <- outer(1:9, 10^(1:5))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1, cex.axis=0.75)
abline(r5, col = "blue", lwd=2)
abline(r5a, col = "olivedrab3", lwd=1, lty=2)
abline(r5p, col = "purple", lwd=1, lty=2)
abline(r5w, col = "orange", lwd=1, lty=2)
mtext("number of cytotypes", side=1, line=2, outer=TRUE, cex=0.55)
mtext("(log10 scale)", side=1, line=3, outer=TRUE, cex=0.55)
mtext("net diversification rate", side=2, line=0, outer=TRUE, cex=0.55, at=0.9)
mtext("net diversification rate", side=2, line=0, outer=TRUE, cex=0.55, at=0.7)
mtext("net diversification rate", side=2, line=0, outer=TRUE, cex=0.55, at=0.5)
mtext("(epsilon=0)", side=2, line=-1, outer=TRUE, cex=0.55, at=0.9)
mtext("(epsilon=0.9)", side=2, line=-1, outer=TRUE, cex=0.55, at=0.7)
mtext("(MEDUSA)", side=2, line=-1, outer=TRUE, cex=0.55, at=0.5)
mtext("extant species richness", side=2, line=0, outer=TRUE, cex=0.55, at=0.3)
mtext("(SR)", side=2, line=-1, outer=TRUE, cex=0.55, at=0.3)
mtext("stem age", side=2, line=0, outer=TRUE, cex=0.55, at=0.1)
mtext("(Myr)", side=2, line=-1, outer=TRUE, cex=0.55, at=0.1)
dev.off()

####
models<-list(r1, r2, r3, r4, r5, r1a, r2a, r3a, r4a, r5a, r1p, r2p, r3p, r4p, r5p, r1w, r2w, r3w, r4w, r5w)
variable<-c("div0", "div0.9", "r.mean", "SR", "age", "div0", "div0.9", "r.mean", "SR", "age", "div0", "div0.9", "r.mean", "SR", "age", "div0", "div0.9", "r.mean", "SR", "age")

model.list = list()
for (i in 1:length(models)){      
  tabellino<-data.frame(variable[i]
                        , paste(round(intervals(models[[i]], which="all")$coef[1,2], 3), " (", round(intervals(models[[i]], which="all")$coef[1,1], 3 ), ", ", round(intervals(models[[i]], which="all")$coef[1,3], 3), ")", sep="")
                        , paste(round(intervals(models[[i]], which="all")$coef[2,2], 3), " (", round(intervals(models[[i]], which="all")$coef[2,1], 3 ), ", ", round(intervals(models[[i]], which="all")$coef[2,3], 3), ")", sep="")
                        , round(summary(models[[i]])$tTable[2,4], 3)
                        , round(R2.pred(mod = models[[i]]), 4)
                        , paste(round(intervals(models[[i]], which="all")$corStruct[1,2], 3), " (", round(intervals(models[[i]], which="all")$corStruct[1,1], 3 ), ", ", round(intervals(models[[i]], which="all")$corStruct[1,3], 3), ")", sep="")
  )
  model.list[[length(model.list)+1]] = tabellino
}

tab.res <- do.call(rbind, model.list)
clades<-c("angiosperms", "angiosperms", "angiosperms", "angiosperms", "angiosperms", "monocots", "monocots", "monocots", "monocots", "monocots", "rosids", "rosids", "rosids", "rosids", "rosids", "asterids", "asterids", "asterids", "asterids", "asterids")
numbers<-c(length(cdat[[1]]$tip.label), length(cdat[[1]]$tip.label), length(cdat[[1]]$tip.label), length(cdat[[1]]$tip.label), length(cdat[[1]]$tip.label), length(ca[[1]]$tip.label), length(ca[[1]]$tip.label), length(ca[[1]]$tip.label), length(ca[[1]]$tip.label), length(ca[[1]]$tip.label), length(cp[[1]]$tip.label), length(cp[[1]]$tip.label),  length(cp[[1]]$tip.label), length(cp[[1]]$tip.label), length(cp[[1]]$tip.label), length(cw[[1]]$tip.label), length(cw[[1]]$tip.label), length(cw[[1]]$tip.label), length(cw[[1]]$tip.label), length(cw[[1]]$tip.label))
tab.res<-cbind(numbers, tab.res)
tab.res<-cbind(clades, tab.res)
names(tab.res)<-c("clade", "n (sampled genera)", "predictor", "intercept", "slope", "P", "Rsq", "lambda")
write.table(tab.res, "results/2.results.tab.genera.SB.csv", row.names=FALSE, sep=",")

