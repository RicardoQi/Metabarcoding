###Github link: https://github.com/RicardoQi/Metabarcoding.git
###This code is written by RicardoQi to analyze the relationship between the invasive species, garlic mustard, and local plant communities.
###It uses Boxplot, RDA, and NMDS to help figure it out.

#Packages installation
install.packages("tidyverse")
install.packages("vegan")
install.packages("lmtest")
install.packages("car")
install.packages("coefplot")

library(tidyverse)
library(vegan)
library(lmtest)
library(car)
library(coefplot)

#Importing the csv file and divide it into two parts, one with all species names, and one with conditions of garlic mustard.
FloristicSurvey <- read.csv("C:\\Users\\qi199\\Documents\\Metabarcoding\\Data\\FloristicSurvey.csv")
MyData<-FloristicSurvey[,1:10]
MyData$Population<-as.factor(MyData$Population)
SpeciesComp<-FloristicSurvey[,11:length(FloristicSurvey)]
rownames(SpeciesComp)<-MyData[,1]

#Two columns created to help regulate.
MyData$Richness<-rowSums(SpeciesComp!=0)
MyData$Shannon<-diversity(SpeciesComp, index = "shannon", MARGIN = 1, base = exp(1))

#mod1 is the overall condition of the garlic mustard.
mod1<-lm(Richness~ Population + Location + Rosettes + Bolting+ Budding + Bud_Flw, data = MyData)
car::qqp(MyData$Richness, "norm")

plot(density(MyData$Richness))

mod1<-mod1.mBud_Flw
summary(mod1)

#Boxplot between richness and location created
ggplot(MyData, aes(y=Richness,x=Location))+
  geom_boxplot(aes(color=Location))+
  labs(y="Richness")

#Convert to matrix for further use
SpeciesComp.hell <- decostand(SpeciesComp, "hell")
var <-MyData[,c(2:6)]
formulaRDA<- rda(SpeciesComp.hell ~ Location+Rosettes+Bolting+Budding+Population, data=var, scale=T)

#Plotting RDA
smry <- summary(formulaRDA)
df1  <- data.frame(smry$sites[,1:2])      
df2  <- data.frame(smry$biplot[,1:2]) 
rda.plot <- ggplot(df1, aes(x=RDA1, y=RDA2)) + 
  geom_text(aes(label=rownames(df1)),size=3,position=position_jitter(width=.2,height=.2)) +
  geom_point(aes(alpha=0.3)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted")  


(formulaRDA.PLOT<-rda.plot +
    geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                 color="red", arrow=arrow(length=unit(0.01,"npc"))) +
    geom_text(data=df2, 
              aes(x=RDA1,y=RDA2,label=rownames(df2),
                  hjust=0.5*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
              color="red", size=4)+ theme(legend.position="none"))

#Writting NMDS codes and plot.
M <- as.matrix(SpeciesComp)
dist_M <- vegdist(M, method = "bray", binary = T)
meta.nmds <- metaMDS(dist_M, k=2, trymax = 1000)
NMDS.data<-data.frame(Location=MyData$Location,Population=MyData$Population) 
NMDS.data$NMDS1<-meta.nmds$points[ ,1] 
NMDS.data$NMDS2<-meta.nmds$points[ ,2] 
ggplot(data = NMDS.data, aes(y = NMDS2, x = NMDS1))+ 
  geom_point( aes(color = Population,shape = Location), size = 1.5,alpha=0.6)
