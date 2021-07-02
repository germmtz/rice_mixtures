####### Germain Montazeaud ######
## Analysis of rice mixture data for the study "Crop mixtures: does niche complementarity hold for belowground resources? an experimental test using rice genotypic pairs"

#######################################
#######################################
## Set Script parameters ##
#######################################
#######################################

trait_file <- getURL("https://raw.githubusercontent.com/germmtz/rice_mixtures/Rice_traits.csv") # input data file 

## Required packages ##
library(agricolae)
library(MASS)
library(lattice)
library(Hmisc)
library(corrplot)
library(RColorBrewer)
library(lsmeans)
library(FactoMineR)
library(stringr)
library(ggplot2)

## Saving R default graphic parameters
dft.par <- par()


#######################################
#######################################
## Data loading ##
#######################################
#######################################


setwd(wd) # Set working directory
riz <- read.table(trait_file, header=T, sep=",") # loading input data file
summary(riz)



#########################################
#########################################
#   DESCRIPTION OF MONOCULTURES  #
#########################################
#########################################

# Selecting only monocultures to characterize genotypes #
rizM <- droplevels(riz[which(riz$Asso=="M"),])
rizM$genophospho <- as.factor(paste(rizM$IDgeno, rizM$Traitement, sep="-"))

# Creating a dataframe for to sum the two biomass in each pot #
rizMpot_biom <- aggregate(cbind(RB_top, BIOM_above, Biovolume, Tillers)~IDpot, data=rizM, FUN = sum) # we need to compute the sum for BIOM_aer and BIOM_R_0_20 as it was measured for each plant
rizMpot_aer <- aggregate(cbind(PH,SLA)~IDpot, data=rizM, FUN = mean) # we need to compute the sum for BIOM_aer and BIOM_R_0_20 as it was measured for each plant

rizMpot <- unique(rizM[,which(colnames(rizM)%in%c("IDpot","IDcouple","Bloc","Treatment","RB_deep"))]) # BIOM_R_20_40 & BIOM_R_40_60 were initially measured for the entire pot, so we just retrieve them
rizMpot <- merge(rizMpot, rizMpot_biom, by="IDpot") ## Both biomass are merged
rizMpot <- merge(rizMpot, rizMpot_aer, by="IDpot") 


# Computing total root biomass ("BIOM_R")
rizMpot$BIOM_R <- rizMpot$RB_top + rizMpot$RB_deep

# Computing Total biomass in the pot ("BIOM")
rizMpot$BIOM <- rizMpot$BIOM_R+rizMpot$BIOM_above

# Computing Root:Shoot ratio ("Root:Shoot")
rizMpot$"Root:Shoot" <- rizMpot$BIOM_R/rizMpot$BIOM_above

# Hybrid variable : combination of genotype identity and phosphorous availability (will be used to study GxP interaction)
rizMpot$genophospho <- as.factor(paste(rizMpot$IDcouple, rizMpot$Treatment, sep="-"))



############################
#### PCA representation of monoculture pot ####
###########################

## Building a data-frame with one value per trait per pot. Certain trait will be concatenated with mean, others with sum
PCA1_tab <- rizMpot[,c("IDpot","Bloc","Treatment", "IDcouple","genophospho", "BIOM_above","BIOM_R", "RB_deep","Root:Shoot", "Biovolume","Tillers")] # previous tab is retrieved

PCA2_tab <- aggregate(cbind(PH, SLA, D_bas, SRL_bas, RTD_bas, RBI_bas, PfR_bas, D_ad, SRL_ad, RTD_ad, RBI_ad, PfR_ad)~IDpot, data=rizM, FUN=mean) # all other traits are summed

# Merging all traits
PCA_tab <- merge(PCA1_tab, PCA2_tab, by="IDpot")


# Creating a new variable identifying each rice subspecies (Japonica or Indica)
PCA_tab$subspecie <- rep(NA, nrow(PCA_tab))

for (i in 1:nrow(PCA_tab)) {
  if (PCA_tab[i,"IDcouple"]=="Pdi-Pdi" | PCA_tab[i,"IDcouple"]=="Ktn-Ktn") {
    PCA_tab[i,"subspecie"] <- "JAP"
  } else {
    PCA_tab[i,"subspecie"] <- "IND"
  }
}

PCA_tab$subspecie <- as.factor(PCA_tab$subspecie)

# Hybrid variable Subspecie X Phosphorous (interaction will be investigate)
PCA_tab$subspeciephospho <- as.factor(paste(PCA_tab$subspecie, PCA_tab$Treatment, sep="-"))

# Aspect columns modification
PCA_tab <- PCA_tab[,c(1:5,(ncol(PCA_tab)-1), ncol(PCA_tab),6:(ncol(PCA_tab)-2))]
row.names(PCA_tab) <- PCA_tab$IDpot
colnames(PCA_tab)[colnames(PCA_tab)%in%c("Treatment","IDcouple","subspecie")] <- c("Phosphorous Level","Genotype","Subspecies")
PCA_tab$Genotype <- as.factor(gsub("-.*","",PCA_tab$Genotype))


### Finishing PCA

## Listing all quantitative variables that will be used for PCA
var_quant <- colnames(PCA_tab)[8:ncol(PCA_tab)][!(colnames(PCA_tab)[8:ncol(PCA_tab)]%in%c("BIOM_above","BIOM_R","Root:Shoot"))]

## Listing all qualitative variables that will be used for PCA
var_qualit <- colnames(PCA_tab)[1:7][!(colnames(PCA_tab)[1:7]%in%c("IDpot"))]
# 
## variable selection
PCA_tab <- PCA_tab[,colnames(PCA_tab)%in%c(var_quant, var_qualit)]

## Running PCA ##
res.PCA <- PCA(PCA_tab, quali.sup=1:length(var_qualit), quanti.sup = 12:16,  graph=F)

## Looking eigen values
barplot(res.PCA$eig[,2], names=paste("Dim",1:nrow(res.PCA$eig)))
res.PCA$eig
## 4 first component explain ~91% of total variance


## Outputting PCA graphs
pdf("All_var_PCA_invis.pdf")

par(mar=c(4.5,5.1,3.1,2.1))
plot(res.PCA, choix="var", axes=1:2, cex=1.57,cex.axis=1.5, las=1, title="",cex.lab=1.5, label="none") # add label="none" to do post editing
plotellipses(res.PCA, keepvar =c(1,2,4), axes=1:2, cex=1.5, level=0, pch=20, xlim=c(-8,8), ylim=c(-8,8))

plot(res.PCA, choix="var", axes=2:3, cex=1.57, cex.axis=1.5, las=1, title="")
plotellipses(res.PCA, axes=2:3, keepvar=c(1,2,4), cex=1.5,level=0, pch=20, xlim=c(-7,6), ylim=c(-6,6))

dev.off()

## Looking for links between qualitative variables (phosphorous, subspecies, ...) and individual PCA coordinates (simple way ANOVA)
dimdesc(res.PCA)


## Building a better representation of individual in the PCA space
PC <- cbind( PCA_tab[,c("Bloc","Phosphorous Level","Genotype")], res.PCA$ind$coord[,c(1,2)])


pdf("Individual_PCA.pdf")
par(mfrow=c(1,1), mar=c(5.1, 5.1, 4.1 ,2.1))
barycentre <- aggregate(cbind(Dim.1, Dim.2)~Genotype, data=PC, FUN = mean)

plot(Dim.2~Dim.1, data=PC, las=1, xlim=c(-5,8), ylim=c(-4,4), type="n", xlab="", ylab="", cex=1.2, cex.axis=1.5, cex.lab=1.5)
abline(h=0, lty=2)
abline(v=0, lty=2)

for (g in levels(PC$Genotype)) {
  segments(barycentre[which(barycentre$Genotype==g),"Dim.1"], barycentre[which(barycentre$Genotype==g),"Dim.2"], PC[which(PC$Genotype==g),"Dim.1"],PC[which(PC$Genotype==g),"Dim.2"], col="grey")
}

points (PC$Dim.1,PC$Dim.2, pch=c(15,16,17,18)[as.numeric(PC$Genotype)], col=c("black","grey76")[as.numeric(as.factor(PC$"Phosphorous Level"))], cex=1.2)


points(barycentre$Dim.1, barycentre$Dim.2, pch=c(15,16,17,18), cex=2, col="gray48")

text(2,1,"Ketan", cex=1.5)
text(1.7,-1.8,"Padi", cex=1.5)
text(-3.2,0.2,"IR64", cex=1.5)
text(-0.5,0.8,"IR64+", cex=1.5)

legend("bottomright",fill=c("grey","black"), border=c("grey76","black"),legend=c("P0","P+"), bty="n", cex=1.5)


dev.off()



############################
############################
## Monoculture trait comparison ##
############################
############################

rizMpot$Genotype <- as.factor(gsub("-.*","",rizMpot$IDcouple))
rizMpot$Treatment <-  factor(rizMpot$Treatment, levels = c("P0", "P+"))

## Hybrid variable Genotype X Phosphorous (used later to sudy GxP interaction)
rizMpot$Genotype_Treatment <- as.factor(paste(rizMpot$Genotype, rizMpot$Treatment, sep="-"))
rizMpot$Genotype_Treatment <- factor(rizMpot$Genotype_Treatment, levels=c("I64-P0","I64-P+","I64+-P0","I64+-P+","Ktn-P0", "Ktn-P+","Pdi-P0","Pdi-P+")) 

var <- c("BIOM_above","PH","Biovolume","Tillers","SLA","RB_top", "RB_deep")
nice_legend <- c("Above-ground Biomass (g)", "Plant height (cm)", expression("Biovolume "*"(m"^3*")"), "Number of tillers", expression("SLA "*"(m².kg"^-1*")"), "Root Biomass 0-20 cm (g)","Root Biomass 20-60 cm (g)")
letters <- c("A)","B)","C)","D)", "E)","F)","G)")


### Outputting trait monocultures traot comparisons:

pdf("Monoculture_traits_comparison.pdf", height=7, width=10)

par(mar=c(3,5,2,1))
par(mfrow=c(2,4))


for (i in 1:length(var)) {


  v <- var[i]
  yaxis <- nice_legend[i]
  letter <- letters[i]
  
  # Statistical model
  mod <- lm(rizMpot[,v] ~ Genotype + Treatment + Genotype_Treatment + Bloc, data=rizMpot)
  anova(mod)
  
  # Multiple comparison (Fisher-LSD)
  test <- LSD.test(mod, "Genotype_Treatment")
  
  
  # Creating a multiple comparison table to store each GxP significan difference
  sigtab <- data.frame(ID=rep(NA, nlevels(rizMpot$Genotype_Treatment)), sig=rep(NA, nlevels(rizMpot$Genotype_Treatment)))
  for (i in 1:nlevels(rizMpot$Genotype_Treatment)) {
    sigtab[i,"ID"] <- levels(rizMpot$Genotype_Treatment)[i]
    sigtab[i,"sig"] <- as.character(test$groups[which(row.names(test$groups)==levels(rizMpot$Genotype_Treatment)[i]),"groups"])
  }
  
  ## Computing mean and se (Standard Error) for each GxP combination
  rizRS_means <- aggregate(rizMpot[,v], by = list(Genotype = rizMpot$Genotype, Treatment=rizMpot$Treatment),FUN = function(x) c(mean = mean(x, na.rm = T), se=sd(x, na.rm=T)/sqrt(length(x))))
  rizRS_means <- do.call(data.frame, rizRS_means)
  
  ## Transform outputted means and se in a proper matrice format (enable to group barplots for P+ and P- , to be called later in barplot() function)
  tabmeans <- tapply(rizRS_means$x.mean, list(rizRS_means$Treatment, rizRS_means$Genotype), function(x) c(x = x))
  
  tabse <- tapply(rizRS_means$x.se, list(rizRS_means$Treatment, rizRS_means$Genotype), function(x) c(x = x))
  
  ## Computing upper limit for y axis
  pltop <- max(rizRS_means$x.mean) + rizRS_means[which(rizRS_means$x.mean==max(rizRS_means$x.mean)),"x.se"]*3
  
  ## generating barplot
  bar <- barplot(height = tabmeans,beside=T, las=1, ylab=yaxis, ylim=c(0,pltop+pltop/6), col=c("grey","black"), cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5, args.legend=list(bty="n","topright",cex=1.5)) # Supplementary argument to draw legend on each plot :, legend.text = T
  
  ## Adding error bars
  segments(bar,tabmeans, bar,tabmeans+tabse, lwd = 2, col="black")
  segments(bar,tabmeans-tabse, bar,tabmeans, lwd = 2, col="white")
  
  ## Adding significant differences with LSD letters
  text(bar, tabmeans+tabse+(pltop/20), sigtab[,"sig"], cex=1.3)
  
  ## Adding letters for figure subsections
  if (v%in%c("PH","Biovolume","BIOM_above","Tillers")) {
  mtext(letter, side=3, at=1,line=-0.5, font=2, cex=1.2)
  } else {
    mtext(letter, side=3, at=1,line=1, font=2, cex=1.2)
    
  }
}

dev.off()



##############################
############## Genotypes' Relative Yields (RY) ############
#############################

## retrieving only polyculture data, as we want to study what happens when genotypes are cultivated with a different neighbor
rizP <- droplevels(riz[which(riz$Asso=="P"),])
rizP$IDgeno <- factor(rizP$IDgeno, levels=c("I64", "I64+","Ktn", "Pdi"))
rizP$IDnei <- factor(rizP$IDnei, levels=c("I64", "I64+","Ktn", "Pdi"))

## Hybrid variables, will be used to study DGE X P and IGE x P interactions 
rizP$genophospho <- as.factor(paste(rizP$IDgeno, rizP$Treatment, sep="-"))
rizP$genophospho <- factor(rizP$genophospho, levels=c("I64-P0","I64-P+","I64+-P0","I64+-P+","Ktn-P0", "Ktn-P+","Pdi-P0","Pdi-P+"))
rizP$neiphospho <- as.factor(paste(rizP$IDnei, rizP$Treatment, sep="-"))
rizP$genoneiphospho <- as.factor(paste(rizP$IDgeno, rizP$IDnei, rizP$Treatment, sep="-"))


#### Characterizing genotypes for their reaction to competition with a stranger : Relative Yield index (RY) ####

### Computing mean monoculture above-ground biomass production (for two plants in the pot) per phosphorous treatment
monoc_aer_means <- aggregate(BIOM_above~Genotype+Treatment, data=rizMpot, FUN=mean)

## Computing RY for each plant 

for (i in 1:nrow(rizP)) {
  geno <- rizP[i,"IDgeno"]
  trt <- rizP[i,"Treatment"]
  rizP[i,"RY"] <- rizP[i,"BIOM_above"]/monoc_aer_means[which(monoc_aer_means$Genotype==geno & monoc_aer_means$Treatment==trt),"BIOM_above"]
}



## Same model as modRY, written with hybrid variables
modRY2 <- lm(RY~IDgeno*IDnei + Treatment + Bloc + genophospho + neiphospho + genoneiphospho, data=rizP)
anova(modRY2)
testRY.perse <- LSD.test(modRY2, "genophospho")
testRY.nei <- LSD.test(modRY2, "neiphospho")

## Tab to store two statics : letters of significance for mutliple comparison tests (testRY.perse, testRY.nei) and significance of a t.test to check if RYs are significantly different from 0.5

RYsig <- data.frame(genophospho=as.character(levels(rizP$genophospho)))
for (c in levels(rizP$genophospho)) {
  RYsig[RYsig$genophospho==c,"mean.perse"] <- mean(rizP[which(rizP$genophospho==c), "RY"], na.rm=T)
  RYsig[RYsig$genophospho==c,"p.val.perse"] <- t.test(rizP[which(rizP$genophospho==c), "RY"], mu=0.5)$p.value
  RYsig[RYsig$genophospho==c,"sig.group.per.se"] <- as.character(testRY.perse$groups[which(row.names(testRY.perse$groups)==c), "groups"])

  RYsig[RYsig$genophospho==c,"mean.nei"] <- mean(rizP[which(rizP$neiphospho==c), "RY"], na.rm=T)
  RYsig[RYsig$genophospho==c,"p.val.nei"] <- t.test(rizP[which(rizP$neiphospho==c), "RY"], mu=0.5)$p.value
  RYsig[RYsig$genophospho==c,"sig.group.nei"] <- as.character(testRY.nei$groups[which(row.names(testRY.nei$groups)==c), "groups"])
  
}

## Generating a graphs with RY profiles 

pdf("RY_profiles.pdf", height=7, width=10)
par(mar=c(4.1,5.1,4.1,6.1))

## PerSe profiles
boxplot(RY~genophospho, data=rizP, ylab="RY", cex.axis=1.5, cex=1.5, cex.lab=1.5, las=1, xaxt="n", at=c(1,2,4,5,7,8,10,11), col=rep(c("grey","black"),4), medcol=rep(c("black","white"),4), ylim=c(0.15,1.25), xlab="")
abline(h=0.5, lty=2)

text(c(1,7,8), c(0.38,0.57,0.57), label="*", cex=3, font=2, col=c("black","black","white"))


## Axis with genotypes Identity
axis(1, at=seq(1.5,10.5, by=3), labels=levels(rizP$IDgeno), cex.axis=1.5)

## Axis with significance for multiple-comparison test (LSD)
axis(3, at=c(1,2,4,5,7,8,10,11), labels=RYsig$sig.group.per.se, cex.axis=1.5,line=-0.16,tcl=0.3, padj=4, font=2, cex=2, lwd=0)

## Legend for P+ and P0
legend("right", xpd=TRUE, inset=c(-0.12,0), legend=c("P0","P+"), fill=c("grey","black"), bty="n", cex=1.5)

dev.off()

###########################################################
## Relationships between RY and hierarchical trait distance (only traits for which the effect of the distance is significant are considered, based on preliminary analysis)
###########################################################

## Computing trait distances
var <- c("Biovolume","SLA","RB_top","D_ad","SRL_ad")
for (v in var) {
  for (i in 1:nrow(rizP)) {
    pot<- rizP[i,"IDpot"]
    plant <- rizP[i,"IDplant"]
    rizP[i,paste("Hdist",v,sep="_")] <- rizP[i,v]-rizP[which(rizP$IDpot==pot & rizP$IDplant!=plant),v]
    }
  }


nice_legend <- c(expression("Biovolume Distance (m"^3*")"),expression("SLA Distance (m"^2*".kg"^-1*")"),expression("RB"["top"]*" Distance (g)"),expression("D"["ad"]*" Distance (mm)"), expression("SRL"["ad"]*" Distance (m.g"^-1*")"))
letters <- c("A)","B)","C)","D)","E)")

var <- paste("Hdist",c("Biovolume","SLA","RB_top","D_ad","SRL_ad"), sep="_")

rizPP0 <- rizP[which(rizP$Treatment=="P0"),]
rizPPp <- rizP[which(rizP$Treatment=="P+"),]

pdf("RY_vs_TraitsDistances.pdf", height=4, width=8, pointsize = 0.5)
par(mfrow=c(2,3),mar=c(4,5,3,3), oma=c(8,8,8,8))
  i <- 1 
  v <-as.character(var[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  mod1 <- lm(RY~rizP[,v]*Treatment+Bloc, data=rizP)
  anova(mod1)
  
  plot(RY~rizPP0[,v], data=rizPP0, las=1, bty="l", ylab="RY", xlab=xlab, pch=16, col=c("grey"), cex.lab=2, cex.axis=2, cex=1.5)
  mtext(letter, side=3, line=0.5,  adj=0, font=2, cex=2)
  points(rizPPp[,v], rizPPp$RY, pch=16, col="black",  cex=1.5)
  
  mod2 <- lm(RY~rizPP0[,v], data=rizPP0)
  mod3 <- lm(RY~rizPPp[,v], data=rizPPp)
  abline(mod2, lty=2, col="grey", lwd=2)
  abline(mod3, lty=2, col="black", lwd=2)
  summary(mod2)
  summary(mod3)
  
  legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2),"", sep=""), paste("R² = ", round(summary(mod3)$r.squared,2),"***", sep="")), bty="n", cex=2, text.col=c("grey","black"))
      
  i <- 2
  v <-as.character(var[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  mod1 <- lm(RY~rizP[,v]*Treatment+Bloc, data=rizP)
  anova(mod1)

  plot(RY~rizPP0[,v], data=rizPP0, las=1, bty="l", ylab="RY", xlab=xlab, pch=16, col=c("grey"), cex.lab=2, cex.axis=2, cex=1.5)
  points(rizPPp[,v], rizPPp$RY, pch=16, col="black",  cex=1.5)
  mtext(letter, side=3, line=0.5, adj=0, font=2, cex=2)
  
  mod2 <- lm(RY~rizP[,v], data=rizP)
  abline(mod2, lty=2, col="black", lwd=2)

  legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=paste("R² = ", round(summary(mod2)$r.squared,2),"***", sep=""), bty="n", cex=2, text.col="black")

  i <- 3
  v <-as.character(var[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  mod1 <- lm(RY~rizP[,v]*Treatment+Bloc, data=rizP)
  anova(mod1)

  plot(RY~rizPP0[,v], data=rizPP0, las=1, bty="l", ylab="RY", xlab=xlab, pch=16, col=c("grey"), cex.lab=2, cex.axis=2, cex=1.5)
  mtext(letter, side=3, line=0.5, adj=0, font=2, cex=2)
  points(rizPPp[,v], rizPPp$RY, pch=16, col="black",  cex=1.5)

  mod2 <- lm(RY~rizPP0[,v], data=rizPP0)
  mod3 <- lm(RY~rizPPp[,v], data=rizPPp)
  abline(mod2, lty=2, col="grey", lwd=2)
  abline(mod3, lty=2, col="black", lwd=2)
  summary(mod2)
  summary(mod3)

  legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2),"*", sep=""), paste("R² = ", sprintf("%.2f",round(summary(mod3)$r.squared,2)),"***", sep="")), bty="n", cex=2, text.col=c("grey","black"))

  i <- 4
  v <-as.character(var[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  mod1 <- lm(RY~rizP[,v]*Treatment+Bloc, data=rizP)
  anova(mod1)

  plot(RY~rizPP0[,v], data=rizPP0, las=1, bty="l", ylab="RY", xlab=xlab, pch=16, col=c("grey"), cex.lab=2, cex.axis=2, cex=1.5)
  mtext(letter, side=3, line=0.5,  adj=0, font=2, cex=2)
  points(rizPPp[,v], rizPPp$RY, pch=16, col="black",  cex=1.5)
  
  mod2 <- lm(RY~rizP[,v], data=rizP)
  abline(mod2, lty=2, col="black")
  summary(mod2)

  legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=paste("R² = ", round(summary(mod2)$r.squared,2),"*", sep=""), bty="n", cex=2, text.col="black")

  i <- 5
  v <-as.character(var[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  mod1 <- lm(RY~rizP[,v]*Treatment+Bloc, data=rizP)
  anova(mod1)

  plot(RY~rizPP0[,v], data=rizPP0, las=1, bty="l", ylab="RY", xlab=xlab, pch=16, col=c("grey"), cex.lab=2, cex.axis=2, cex=1.5)
  mtext(letter, side=3, line=0.5,  adj=0, font=2, cex=2)
  points(rizPPp[,v], rizPPp$RY, pch=16, col="black",  cex=1.5)
  mod2 <- lm(RY~rizPP0[,v], data=rizPP0)
  abline(mod2, lty=2, col="black")
  summary(mod2)


  legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=paste("R² = ", round(summary(mod2)$r.squared,2),"*", sep=""), bty="n", cex=2, text.col="black")
  
  
dev.off()



#########################################
#########################################
#   Comparing mixtures and monoculture biomass #
#########################################
#########################################

## Building a tab for pot biomass production, either monocultures and polycultures
riz_pot_tot_1 <- unique(riz[,colnames(riz)%in%c("IDpot","IDcouple","Asso","Bloc","Treatment","RB_deep")]) 
riz_pot_tot_2 <- aggregate(cbind(BIOM_above, RB_top)~IDpot, data=riz, FUN = sum)
riz_pot_tot <- merge(riz_pot_tot_1, riz_pot_tot_2, by="IDpot")
riz_pot_tot$BIOM <- riz_pot_tot$BIOM_above + riz_pot_tot$RB_top + riz_pot_tot$RB_deep
riz_pot_tot$Treatment <- factor(riz_pot_tot$Treatment, levels=c("P0", "P+"))

# Hybrid variable : GenoCouple X Phosphorous
riz_pot_tot$genophospho <- as.factor(paste(riz_pot_tot$IDcouple, riz_pot_tot$Treatment, sep="-"))

## Testing if the different mixtures have different productivity (in interaction with P)
mod_mixture_tot <- lm(BIOM~IDcouple+Treatment+Bloc+genophospho, data=riz_pot_tot)
anova(mod_mixture_tot)
summary(mod_mixture_tot)

## Multiple comparisons test
test_mixture_tot <- LSD.test(mod_mixture_tot, "genophospho")
test_mixture_tot_2 <- LSD.test(mod_mixture_tot, "IDcouple")
mod_test <- lm(BIOM~Asso*Treatment+Bloc, data=riz_pot_tot)
anova(mod_test)

## Computing means and se (Standard Error) for each mixture x P combination
rizMX_means_tot <- aggregate(riz_pot_tot$BIOM, by = list(Mixture = riz_pot_tot$IDcouple, Treatment=riz_pot_tot$Treatment), FUN = function(x) c(mean = mean(x, na.rm = T), se=sd(x, na.rm=T)/sqrt(length(x))))

rizMX_means_tot <- do.call(data.frame, rizMX_means_tot)
colnames(rizMX_means_tot) <- gsub("^x", "BIOM",colnames(rizMX_means_tot))

## Adding significant differences between each mixture x P combination
for (i in 1:nrow(rizMX_means_tot)) {
  Mix <- rizMX_means_tot[i,"Mixture"]
  trt <- rizMX_means_tot[i,"Treatment"]
  rizMX_means_tot[i, "sig"] <- test_mixture_tot$groups[which(row.names(test_mixture_tot$groups)==paste(Mix,trt,sep="-")),"groups"]
}

## Outputting the total  biomass barplots
pdf("Total_biomass.pdf", height=4, width = 8)

par(mar=c(4,4.5,2,1))

## Giving the right matrix format to means and se, so that barplot() function will group P+ and P0 biomass for each mixture
tabmeans <- tapply(rizMX_means_tot$BIOM.mean, list(rizMX_means_tot$Treatment, rizMX_means_tot$Mixture), function(x) c(x = x))

tabse <- tapply(rizMX_means_tot$BIOM.se, list(rizMX_means_tot$Treatment, rizMX_means_tot$Mixture), function(x) c(x = x))

## Means are sorted decreasingly (on the basis of a P+ & P- mean)
tabmeans <- as.data.frame(tabmeans)
avg <- apply(tabmeans, 2, mean)
tabmeans <- rbind(tabmeans, Avg=avg)
tabmeans <- tabmeans[,order(tabmeans[nrow(tabmeans),], decreasing = T)]
tabmeans <- as.matrix(tabmeans[-3,])

## An "artefact" table is generated to superpose a specific plot for monocultures, enabling to diffrenciate them from the others
tabmeans_monoc <- tabmeans
tabmeans_monoc[,c(2:4,7:9)] <- NA

## Computing upper limit for y axis
pltop <- max(rizMX_means_tot$BIOM.mean) + rizMX_means_tot[which(rizMX_means_tot$BIOM.mean==max(rizMX_means_tot$BIOM.mean)),"BIOM.se"]*3

## Generating barplot 

bar <- barplot(height = tabmeans, beside=T ,las=1, ylab="Total Biomass (g)", ylim=c(0,pltop), col=c("grey","black"), axisnames=F, cex.axis = 1.2, cex.lab=1.2, names.arg = F)


## Superposing specific barplot for monocultures (hatched)
barplot(height=tabmeans_monoc, add=T, col=c("black","white"), density=20, angle=45, beside=T, ylim=c(0,pltop), lwd=45, axes=F, axisnames=F)

## Adding error bars
segments(bar,tabmeans, bar,tabmeans+tabse, lwd = 2, col="black")
segments(bar,tabmeans-tabse, bar,tabmeans, lwd = 2, col="white")


## Adding legend for hatched monocultures plot
legend(24.5,26,bty="n", density=20, angle=45, fill="gray48", legend="Monocultures")
legend(28,21,bty="n", fill=c("grey","black"), legend=c("P0","P+"))
text(2.5, 26, label="*", cex=1.5, font=2)
## Axis with genotypic mixtures identity
axis(1, at=seq(from=2, to=2+9*3, by=3), tick=F, labels = row.names(test_mixture_tot_2$groups), cex.axis=0.75, line=-0.5)
## Axis with multiple comparison signficance letters (LSD test)
axis(1, at=seq(from=2, to=2+9*3, by=3), labels = test_mixture_tot_2$groups$groups, font=2, cex=2, tick=F, line=1)

dev.off()



##############################
#### COMPUTING RYT : comparison of total mixture biomass production with the monoculture performance : RYT > 1 means overyielding
##############################

## Building a table to store RYT values
rizPpot <- unique(rizP[,colnames(rizP)%in%c("IDpot","Bloc","Asso","Treatment","IDcouple")])

rizPpot_RYT <- aggregate(rizP$RY~rizP$IDpot, data=rizP, FUN = sum) # RYT is computed as the sum of the two RY in each pot
colnames(rizPpot_RYT) <- c("IDpot","RYT")
rizPpot <- merge(rizPpot, rizPpot_RYT, by="IDpot")

## Re-ordering mixture couple levels according to their rank in total biomass production
rizPpot$IDcouple <- factor(rizPpot$IDcouple, levels= c("I64+-Ktn","I64-Ktn","Ktn-Pdi","I64-I64+","I64+-Pdi","I64-Pdi"))

## Hybrid G x P variable
rizPpot$genophospho <- as.factor(paste(rizPpot$IDcouple, rizPpot$Treatment, sep="-"))
rizPpot$genophospho <- factor(rizPpot$genophospho, levels= c("I64+-Ktn-P0","I64+-Ktn-P+","I64-Ktn-P0","I64-Ktn-P+","Ktn-Pdi-P0","Ktn-Pdi-P+" ,"I64-I64+-P0","I64-I64+-P+","I64+-Pdi-P0","I64+-Pdi-P+","I64-Pdi-P0","I64-Pdi-P+"))
rizPpot$Treatment <- factor(rizPpot$Treatment, levels=c("P0", "P+"))

## Look if we have significant differences between mixtures for their RYT, and interactions with Phosphorous
mod_mixture_RYT <- lm(RYT~IDcouple+Treatment+genophospho+Bloc, data=rizPpot)
anova(mod_mixture_RYT)


## Multiple comparison test (Fisher - LSD)
testRYT <- LSD.test(mod_mixture_RYT, "genophospho")

## Building a table to store two different stats : multiple comparison, and t.test() to test if RYT is significantly different from 1
RYTsig <- data.frame(genophospho=as.character(levels(rizPpot$genophospho)))
for (c in levels(rizPpot$genophospho)) {
  RYTsig[RYTsig$genophospho==c,"mean"] <- mean(rizPpot[which(rizPpot$genophospho==c), "RYT"], na.rm=T)
  RYTsig[RYTsig$genophospho==c,"p.val"] <- t.test(rizPpot[which(rizPpot$genophospho==c), "RYT"], mu=1)$p.value
  RYTsig[RYTsig$genophospho==c,"sig.group"] <- as.character(testRYT$groups[which(row.names(testRYT$groups)==c), "groups"])
}


## Setting default graphic paramters
par(dft.par)

## Outputting RYT profiles
pdf("RYT_profiles.pdf", height=7, width=10)

par(mar=c(4.1,5.1,4.1,6.1))

boxplot(RYT~genophospho, data=rizPpot, ylab="RYT", cex.axis=1.5, cex=1.5, cex.lab=1.5, las=1, xaxt="n", at=c(1,2,4,5,7,8,10,11, 13,14, 16,17), col=rep(c("grey","black"),6), medcol=rep(c("black","white"),6), xlab="")
abline(h=1, lty=2)

## Axis with mixtures identify
axis(1, at=seq(1.5,16.5, by=3), labels=levels(rizPpot$IDcouple), cex.axis=1.5)

## Axis with multiple comparison letters
axis(3, at=c(1,2,4,5,7,8,10,11, 13,14, 16,17), labels=RYTsig$sig.group, cex.axis=0.9, line=-0.25, tcl=0.3, padj=4, font=2, lwd=0)

## Legend with phosphorous treatment
legend("right", xpd=TRUE, inset=c(-0.12,0), legend=levels(rizPpot$Treatment), fill=c("grey","black"), bty="n", cex=1.5)


dev.off()


##############################
#### COMPUTING Dmax : comparison of total mixture biomass production with the best component production (in monoculture)
##############################
### Here, total biomass will be used
var <- "BIOM"

## Computing mean monocultures pruduction per pot, per mixture and per treatment
monocmeans <- aggregate(riz_pot_tot[which(riz_pot_tot$Asso=="M"),var]~IDcouple+Treatment, FUN=mean, data=riz_pot_tot[which(riz_pot_tot$Asso=="M"),])
monocmeans$IDgeno <- gsub(".*-","", monocmeans$IDcouple)
colnames(monocmeans)[3] <- var
levels(monocmeans$Treatment) <- c("P0","P+")

## Preparing a tab to store dmax values
Dmaxtab <- riz_pot_tot[which(riz_pot_tot$Asso=="P"),c("IDpot",var)] # retrieving biomass produced by the two plants in polyculture condition
Dmaxtab <- droplevels(unique(merge(Dmaxtab, rizPpot[,c("IDpot","IDcouple","Treatment","Bloc")], by="IDpot"))) # reporting index columns


# loop to compute Dmax, as the difference between mixtures total biomass and the biomass of the most productive component of the mixture #

for (i in 1:nrow(Dmaxtab)) {
  trt <- Dmaxtab[i, "Treatment"]
  geno <- gsub("-.*", "", Dmaxtab[i,"IDcouple"])
  nei <- gsub(".*-", "", Dmaxtab[i,"IDcouple"])
  Dmaxtab[i,"Dmax"] <- (Dmaxtab[i,var]-max(monocmeans[which(monocmeans$Treatment==trt & gsub(".*-","",monocmeans$IDgeno)%in%c(geno, nei)),var]))/max(monocmeans[which(monocmeans$Treatment==trt & gsub(".*-","",monocmeans$IDgeno)%in%c(geno, nei)),var])
  
}

## Checking variables explaining Dmax
modDmax <- lm(Dmax~IDcouple*Treatment+Bloc, data=Dmaxtab)
anova(modDmax)

## Hybrid Mixture X P variable, to investigate M x P 
Dmaxtab$Couplephospho <- as.factor(paste(Dmaxtab$IDcouple, Dmaxtab$Treatment, sep="-"))

## Reordering Mixture x P variable levels
Dmaxtab$Couplephospho <- factor(Dmaxtab$Couplephospho, levels=c("I64+-Ktn-P0","I64+-Ktn-P+", "I64-Ktn-P0", "I64-Ktn-P+", "Ktn-Pdi-P0", "Ktn-Pdi-P+", "I64-I64+-P0", "I64-I64+-P+", "I64+-Pdi-P0", "I64+-Pdi-P+","I64-Pdi-P0", "I64-Pdi-P+"))

## Same model as modDax, written with hybrid variable
modDmax <- lm(Dmax~IDcouple+Treatment+Bloc+Couplephospho, data=Dmaxtab)
## Multiple comparison (Fisher-LSD) test
testDmax <- LSD.test(modDmax, "Couplephospho")


## Table to store two different statistics : multiple comparison of Dmax for each Mixture x P combination and t.test() to check if Dmax is significantly different from 0 for each mixture x p combination

Dmaxsig <- data.frame(Couplephospho=as.character(levels(Dmaxtab$Couplephospho)))
for (c in levels(Dmaxtab$Couplephospho)) {
  Dmaxsig[Dmaxsig$Couplephospho==c,"mean"] <- mean(Dmaxtab[which(Dmaxtab$Couplephospho==c), "Dmax"], na.rm=T)
  Dmaxsig[Dmaxsig$Couplephospho==c,"p.val"] <- t.test(Dmaxtab[which(Dmaxtab$Couplephospho==c), "Dmax"])$p.value
  Dmaxsig[Dmaxsig$Couplephospho==c,"sig.group"] <- testDmax$groups[which(row.names(testDmax$groups)==c),"groups"]

}


## Outputting Dmax profiles ##

pdf("Dmax_profiles.pdf", height=7, width=10)

par(mar=c(4.1,5.1,4.1,6.1))

boxplot(Dmax~Couplephospho, data=Dmaxtab, ylab="Dmax", cex.axis=1.5, cex=1.5, cex.lab=1.5, las=1, xaxt="n", at=c(1,2,4,5,7,8,10,11,13,14,16,17), col=rep(c("grey","black"),6), medcol=rep(c("black","white"),6), ylim=c(-0.65,0.45), xlab="")
abline(h=0, lty=2)

text(c(4,7,11,14,16), c(-0.2,-0.2,-0.15,-0.27,-0.3), label="*", cex=2.5, font=2, col=c("black","black","white", "white","black"))

## Axis with mixtures identify
axis(1, at=seq(1.5,16.5, by=3), labels=levels(Dmaxtab$IDcouple), cex.axis=1.5)

## Axis with mutliple comparison letters (LSD test)
axis(3, at=c(1,2,4,5,7,8,10,11, 13,14, 16,17), labels=Dmaxsig$sig.group, cex.axis=0.9, line=-0.25,tcl=0.3, padj=4, font=2, lwd=0)

## Legend for phosphorous treatement
legend("right", xpd=TRUE, inset=c(-0.12,0), legend=levels(Dmaxtab$Treatment), fill=c("grey","black"), bty="n", cex=1.5)

dev.off()



#############################
#############################
### Relationships between relative mixture productivity and absolute phenotypic distances
#############################
#############################
rizpt <- riz_pot_tot[,c(1:4,7,10)]
rizpt <- droplevels(rizpt[rizpt$Asso=="P",])

var <- c("Tillers", "PH","Biovolume","SLA","RB_top","D_bas","SRL_bas","RTD_bas","RBI_bas","PfR_bas","D_ad","SRL_ad","RTD_ad","RBI_ad", "PfR_ad", "RY") ## choosing variables for which we will calculate phenotypic distances

## Loop to calculate absolute distances on the selected traits

for (v in var) {
  for (i in 1:nrow(rizpt)) {
    pot <- rizpt[i,"IDpot"]
    rizpt[i,paste("Adist", v, sep="_")] <- abs(rizP[which(rizP$IDpot==pot),v][1]-rizP[which(rizP$IDpot==pot),v][2])
    
  }
}


summary(rizpt)

## Statistical testing of the relation between genotype dissimilarity (absolute phenotypic distance) and mixture Dmax and RYT

## Building the right tab with both Dmax and RYT indexes
riz_mix_index <- merge(rizPpot[,-c(7)], Dmaxtab[,c("IDpot","Dmax")], by="IDpot")
riz_mix_index <- merge(riz_mix_index, rizpt[,grep("IDpot|Adist", colnames(rizpt))], by="IDpot")


## Statistical testing of the relation between genotype dissimilarity and RYT
###############################################################################

## Empty tab to store statistics (degree of freedom, F value, p-value, etc)
stattab <- data.frame(Effect=NULL, Df= NULL, SumSq = NULL, MeanSq=NULL, Fvalue=NULL, Pvalue=NULL)
## Creating an empty row that will seperate each trait distance statistics in the final inegrative table
separation_row <- data.frame(Effect=NA, Df= NA, SumSq = NA, MeanSq=NA, Fvalue=NA, Pvalue=NA)

## Loop to check relation between RYT and phenotypic distance for each trait
for (v in grep("Adist", colnames(riz_mix_index), value=T)) {
  mod <- lm(RYT~riz_mix_index[,v]*Treatment+Bloc, data=riz_mix_index)
  stat <- as.data.frame(anova(mod))
  colnames(stat) <- c("Df","SumSq","MeanSq","Fvalue","Pvalue")
  row.names(stat) <- NULL
  stat <- cbind(Effect=c(v, "Treatment","Bloc", paste(v,":Treatment", sep=""), "Residuals"), stat)
  ## printing each significant effect
  print(stat[stat$Pvalue<0.05,c("Effect","Pvalue")][!is.na(stat[stat$Pvalue<0.05,c("Effect","Pvalue")])])
  
  stat <- rbind(stat, separation_row)
  
  stattab <- rbind(stattab, stat)
  
}
## Outputting statistics in a single .csv file
write.csv(stattab, file="RYT_vs_distances_Stats.csv",row.names=F)

## Exatracting significant effect to generate plots only for effects that have a relation with RYT
significant <- stattab[stattab$Pvalue<0.05, "Effect"]
significant <- significant[!is.na(significant)][-grep("\\:",significant[!is.na(significant)])]
significant <- significant[grep("Adist",significant)]


## Outputting interesting relations between RYT and phenotypic distances

nice_legend <- c(expression("Biovolume Distance (m"^3*")"),expression("SLA Distance (m"^2*".kg"^-1*")"),expression("RB"["top"]*" Distance (g)"),expression("D"["ad"]*" Distance (mm)"), expression("SRL"["ad"]*" Distance (m.g"^-1*")"), "Relative Yield (RY) Distance")
letters <- c("A)","B)","C)","D)","E)","F)")


pdf("RYT_vs_PhenoDist.pdf", height=4, width=8, pointsize = 0.5)
par(mfrow=c(2,3),mar=c(4,5,3,3), oma=c(8,8,8,8))

for (i in 1:length(significant)) {
  
  v <- as.character(significant[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  
  plot(RYT~riz_mix_index[,v], data=riz_mix_index, las=1, bty="l",xlab=xlab,ylab="",pch=16, col=c("grey","black")[as.numeric(riz_mix_index$Treatment)],cex.axis=2, cex.lab=2, cex=1.5, axes=F)
  
  axis(1, cex.lab=2, cex.axis=2)
  axis(2, las=1, cex.lab=2, cex.axis=2)
  mtext(text=expression ("RYT"), side=2, line=4, cex=2)
  box(bty="l")
  
  mod2 <- lm(RYT~riz_mix_index[,v], data=riz_mix_index)
  abline(mod2, lty=2, lwd=2)
  ## Legend with R² and p-value for the two regressions
  
  mtext(letter, side=3, line=0.5, at=0, adj=0.5, font=2, cex=2)
  
  if (summary(mod2)$coef[2,1]>0) {
    
    if (anova(mod2)$"Pr(>F)"[1]<0.001) {
      
      legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),"p <0.001"), bty="n", cex=2)
      
    } else {
      
      legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),paste("p = ", round(anova(mod2)$"Pr(>F)"[1],3))), bty="n", cex=2)
      
    }
    
  } else {
    
    if (anova(mod2)$"Pr(>F)"[1]<0.001) {
      
      legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),"p <0.001"), bty="n", cex=2)
      
    } else {
      
      legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),paste("p = ", round(anova(mod2)$"Pr(>F)"[1],3))), bty="n", cex=2)
      
    }
    
  }
  print(v)
}
dev.off()



## Statistical testing of the relation between phenotypic distances and Dmax
###############################################################################
## Empty tab to store statistics (degree of freedom, F value, p-value, etc)
stattab <- data.frame(Effect=NULL, Df= NULL, SumSq = NULL, MeanSq=NULL, Fvalue=NULL, Pvalue=NULL, slope=NULL)
## Creating an empty row that will seperate each trait distance statistics in the final inegrative table
separation_row <- data.frame(Effect=NA, Df= NA, SumSq = NA, MeanSq=NA, Fvalue=NA, Pvalue=NA, slope=NA)

## Loop to check relation between Dmax and phenotypic distance for each trait
for (v in grep("Adist", colnames(riz_mix_index), value=T)) {
  ## Statistical model
  mod <- lm(Dmax~riz_mix_index[,v]*Treatment+Bloc, data=riz_mix_index)
  stat <- as.data.frame(anova(mod))
  colnames(stat) <- c("Df","SumSq","MeanSq","Fvalue","Pvalue")
  row.names(stat) <- NULL
  stat <- cbind(Effect=c(v, "Treatment","Bloc", paste(v,":Treatment", sep=""), "Residuals"), stat, slope=rep(NA,5))
  stat[1,"slope"] <- summary(mod)$coef[2,1]
  ## printing each significant effect
  print(stat[stat$Pvalue<0.05,c("Effect","Pvalue")][!is.na(stat[stat$Pvalue<0.05,c("Effect","Pvalue")])])
  
  stat <- rbind(stat, separation_row)
  
  stattab <- rbind(stattab, stat)
  
}
## Outputting statistics in a single .csv file
write.csv(stattab, file="Dmax_vs_distances_Stats.csv",row.names=F)

## Exatracting significant effect to generate plots only for effects that have a relation with Dmax
significant <- stattab[stattab$Pvalue<0.05, "Effect"]
significant <- significant[!is.na(significant)]
significant <- significant[grep("Adist",significant)]

## Ouputting interesting relations between Dmax and phenotypic dissimilarities


nice_legend <- c(expression("Biovolume Distance (m"^3*")"),expression("SLA Distance (m"^2*".kg"^-1*")"),expression("RB"["top"]*" Distance (g)"),expression("D"["ad"]*" Distance (mm)"), expression("SRL"["ad"]*" Distance (m.g"^-1*")"), "Relative Yield (RY) Distance")
letters <- c("A)","B)","C)","D)","E)","F)")

pdf("Dmax_vs_PhenoDist.pdf", height=4, width=8, pointsize = 0.5)

par(mfrow=c(2,3), mar=c(4,5,3,3), oma=c(8,8,8,8))

for (i in c(1:(length(significant)))) {
  v <- as.character(significant[i])
  xlab <- nice_legend[i]
  letter <- letters[i]
  
  plot(Dmax~riz_mix_index[,v], data=riz_mix_index, las=1, bty="l",xlab=xlab,ylab="",pch=16, col=c("grey","black")[as.numeric(riz_mix_index$Treatment)],cex.axis=2, cex.lab=2, cex=1.5, axes=F)
  
  axis(1, cex.lab=2, cex.axis=2)
  axis(2, las=1, cex.lab=2, cex.axis=2)
  mtext(text=expression ("D"["max"]), side=2, line=4, cex=2)
  box(bty="l")
  
  mod2 <- lm(Dmax~riz_mix_index[,v], data=riz_mix_index)
  abline(mod2, lty=2, lwd=2)
  ## Legend with R² and p-value for the two regressions
  
  mtext(letter, side=3, line=0.5, at=0, adj=0.5, font=2, cex=2)
  
  if (summary(mod2)$coef[2,1]>0) {
    
    if (anova(mod2)$"Pr(>F)"[1]<0.001) {
      
      legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),"p <0.001"), bty="n", cex=2)
      
    } else {
      
      legend("bottomright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),paste("p = ", round(anova(mod2)$"Pr(>F)"[1],3))), bty="n", cex=2)
      
    }
    
  } else {
    
    if (anova(mod2)$"Pr(>F)"[1]<0.001) {
      
      legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),"p <0.001"), bty="n", cex=2)
      
    } else {
      
      legend("topright", xpd=TRUE, inset=c(-0.05,0), legend=c(paste("R² = ", round(summary(mod2)$r.squared,2)),paste("p = ", round(anova(mod2)$"Pr(>F)"[1],3))), bty="n", cex=2)
      
    }
    
  }
  print(v)
}

dev.off()

