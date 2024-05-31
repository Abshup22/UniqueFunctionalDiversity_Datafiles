library(ggrepel)
library(FD)
library(tidyr)  #nEEDED FOR UNITE
library(ggplot2)
library(ISLR)
library(tree)
library(rpart)
library(rpart.plot)
library(ggpmisc)
library(dplyr)
library(fuzzySim)#USED FOR MAKING PRESENCE-ABSENCE DATA
library(reshape)
library(devtools)
library(vegan)
library(plyr)

#START
####################################################################################
#Uploading database, the database modified from Kate that has added communities
#and all traits
PaleoShup<-read.csv("NAMammalTraits_Database.csv", header = TRUE)

#To compute FD indices, you must have a species-by-species trait matrix or distance matrix

#Only using Name and traits, use columns (LogMass, Locomotion, Life Habit and Recode.Diet)
PPShup1<-as.data.frame(PaleoShup[ ,c(11,13,15, 16, 18)])

duplicated(PPShup1)
PPShup1[duplicated(PPShup1), ]  #Lists all of the repeated data
Paleo1<-PPShup1[!duplicated(PPShup1), ]    #Removes duplicated data and lists one of each, (1759 rows)                

#Named the species as the row names for the traits
rownames(Paleo1)<-Paleo1$Name #Cannot have duplicate names, need to get rid of duplicates


Paleo1 <- Paleo1[order(rownames(Paleo1)),] #making rownames the species name
Paleoy<- Paleo1[,-1]
#Putting rows in alphabetical order to match abun

#################################################################################
#Creating presence-absence data matrix, localities are sites, name is sp.col
PaleoShup<-read.csv("NAMammalTraits_Database.csv", header = TRUE, check.names = F)

#Must create a presence-absence table using the sites by species for dbFD
Paleoabun<-as.data.frame(PaleoShup[ ,c(2,11)]) #isolating species and sites, sites=y, species=x

#Creating a presence-absence data base creating off of data
abun<-splist2presabs(Paleoabun, sites.col="Locality", sp.col="Name", keep.n = TRUE)
rownames(abun)<-abun$Locality
abun1<-abun[,-1]
matabun<-as.matrix(abun1) #needs to be a matrix

################################################################################
#Next step is to run the matrix and dataframe together 

#Finding the rownames and colnames that don't match
row.names(Paleoy)[!(row.names(Paleoy) %in% row.names(t(matabun)))]
row.names(t(matabun))[!(row.names(t(matabun)) %in% row.names(Paleoy))]  #allowing me to find the culprit of non-unqiue

#checking the rownames to compare (FALSE)
identical(row.names(Paleoy), row.names(t(matabun))) 

#Separating the names not known and combining them so that we can figure out the differences
aa <- (1:nrow(Paleoy))[!rownames(Paleoy) %in% colnames(abun1)]
zz <- (1:ncol(abun1))[!colnames(abun1) %in% rownames(Paleoy)]
cbind(rownames(Paleoy)[aa],colnames(abun1)[zz])


colnames(matabun)<-rownames(Paleoy) 
rownames(Paleoy)<-colnames(matabun) #This made the rownames and the colnames the same
#of the data.frame and matrix

#checking the rownames to compare (TRUE)
identical(row.names(Paleoy), row.names(t(matabun))) 

###########################################################################################################
###PCOA 
gow<-gowdis(Paleoy)########Make distance matrix

#####################################################################


ex2<-dbFD(gow, matabun,  corr = "sqrt", stand.x=FALSE, stand.FRic = TRUE, w.abun = FALSE, calc.FRic = TRUE, 
          calc.FDiv = TRUE, m=5, print.pco = TRUE)    

################################################################################################################
#SAVING INFORMATION
write.csv(ex2$FRic, file = "FRIC.csv")
write.csv(ex2$FEve, file = "FEVE.csv")
write.csv(ex2$FDiv, file = "FDIV.csv")
write.csv(ex2$x.values, file="eigenvalues.csv")
write.csv(ex2$x.axes, file = "Speciescoords.csv")

#######################################################################################################

#GETTING MEAN and STDEV of TIME BINS.  1 million year time bins
#Then will compare to the mean and stdev. found above for each time bin
Zachos<-read.csv("Westerhold et al. 2020 (Zachos Data).csv")

Wagner<-read.csv("Richness_per_Myr_Bin_Wagner.csv")
Time.Bins<-as.data.frame(matrix(nrow=66,ncol=17))

colnames(Time.Bins)<-c("time_bin","Mean","SD", "Prop.Archaic")
Time.Bins$time_bin<-1:66

All<-read.csv("NAMammalTraits_Database.csv")

##### Loop start ####
for(i in Wagner$Ma){ 
  #### Data Formatting #### 
  Bin.Name<-i

c<- Zachos[Zachos$age_tuned>= i-1 & Zachos$age_tuned<= i , ] 
(1:length(c$ISOBENbinned_d18O))[is.na(c$ISOBENbinned_d18O)]    #Are there any NA's
Average<-mean(as.numeric(c$ISOBENbinned_d18O[!is.na(c$ISOBENbinned_d18O)]))
SD<-sd(c$ISOBENbinned_d18O[!is.na(c$ISOBENbinned_d18O)])

b<- All[All$Ma>= i-1 & All$Ma<= i , ] 

Prop<-(colSums(b==0)/nrow(b))*100

Prop<-as.data.frame(Prop)
Prop.Archaic<-Prop[13,1]

#fill the dataframe
Time.Bins[i,1] <- Bin.Name # name
Time.Bins[i,2] <- Average # name
Time.Bins[i,3] <- SD # name
Time.Bins[i,4] <- Prop.Archaic

}

#calculatin average FD of each metric for each time bin
Communities<-read.csv("FunctionalCommunities_Final.csv") 
Feve<- Communities%>%
 dplyr:: group_by(Communities$Mil_TimeBin)%>%
  dplyr::summarize(Mean = mean(Feve, na.rm=TRUE))


Feve<-as.data.frame(Feve)
Feve$Ma<-Feve$`Communities$Mil_TimeBin`
Feve$Feve_mean<-Feve$Mean

Communities<-read.csv("FunctionalCommunities_Final.csv") #Find average prop.of archaic mammals
Local_Prop<- Communities%>%
  dplyr:: group_by(Communities$Mil_TimeBin)%>%
  dplyr::summarize(Mean = mean(Prop.Archaic, na.rm=TRUE))
Local_Prop<-as.data.frame(Local_Prop)
write.csv(Local_Prop, file='Local_Prop.csv')


#combining data for 1 myr time bins

Wagner1<-full_join(Wagner1, Feve, by = "Ma")

Wagner$do18<-Time.Bins$Mean[match(Time.Bins$time_bin, Wagner$Ma)]

Wagner$do18_sd<-Time.Bins$SD[match(Time.Bins$time_bin, Wagner$Ma)]

Wagner$Prop.Archaic<-Time.Bins$Prop.Archaic[match(Time.Bins$time_bin, Wagner$Ma)]

write.csv(Wagner1, file = '1myr_TimeBin_Data_Mixed_Model.csv')


######################################################################################

#Mixed Model
MixedModel_Data<-read.csv("1myr_TimeBin_Data_Mixed_Model.csv")

library(lme4)
#Using average local FD indices
lm1<-glm(Fric_mean~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)
lm1$coefficients
plot(lm2)

lm2<-glm(Fric_mean~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2)

lm3<-glm(Fric_mean~MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)

lm4<-glm(Fric_mean~+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Fric_mean~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm6<-glm(Fric_mean~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data) #species richness lognormal
summary(lm1)

lm7<-glm(Fric_mean~MixedModel_Data$lognormal+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #species richness lognormal
summary(lm1)
#########################################################33
#Fdiv
lm1<-glm(Fdiv_mean~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)
plot(lm2)

lm2<-glm(Fdiv_mean~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2)

lm3<-glm(Fdiv_mean~MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)

lm4<-glm(Fdiv_mean~+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Fdiv_mean~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm1<-glm(Fdiv_mean~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data) #species richness lognormal
summary(lm1)

lm1<-glm(Fdiv_mean~MixedModel_Data$lognormal+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #species richness lognormal
summary(lm1)


#######################################################################
#Feve 
lm1<-glm(Feve_mean~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)
plot(lm2)

lm2<-glm(Feve_mean~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2)

lm3<-glm(Feve_mean~MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)

lm4<-glm(Feve_mean~+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Feve_mean~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm1<-glm(Feve_mean~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data) #species richness lognormal
summary(lm1)

lm1<-glm(Feve_mean~MixedModel_Data$lognormal+MixedModel_Data$Local_PropArchaic_Mean, data=MixedModel_Data) #species richness lognormal
summary(lm1)





#Using FD indices for the continental fauna -____________________Fric
library(lme4)
lm1<-glm(Faunal_Fric~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)
plot(lm1)

lm2<-glm(Faunal_Fric~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2) ##****************
plot(lm2)

lm3<-glm(Faunal_Fric~MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)  #********************** Smallest AIC

lm4<-glm(Faunal_Fric~+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Faunal_Fric~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm5<-glm(Faunal_Fric~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm1<-glm(Faunal_Fric~MixedModel_Data$lognormal+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #species richness lognormal
summary(lm1)

#Using FD indices for the continental fauna - -------- Fdiv
library(lme4)
lm1<-glm(Faunal_Fdiv~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)#***************
plot(lm1)

lm2<-glm(Faunal_Fdiv~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2) ##***************************** Smallest AIC
plot(lm2)

lm3<-glm(Faunal_Fdiv~MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)  #*********

lm4<-glm(Faunal_Fdiv~+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Faunal_Fdiv~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm5<-glm(Faunal_Fdiv~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm1<-glm(Faunal_Fdiv~MixedModel_Data$lognormal+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #species richness lognormal
summary(lm1)




#Using FD indices for the continental fauna - -------- Feve
library(lme4)
lm1<-glm(Faunal_Feve~MixedModel_Data$lognormal, data=MixedModel_Data) #species richness lognormal
summary(lm1)#***************
plot(lm1)

lm2<-glm(Faunal_Feve~MixedModel_Data$do18, data=MixedModel_Data) #d018
summary(lm2) ##***************************** Smallest AIC
plot(lm2)

lm3<-glm(Faunal_Feve~MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #proportion of archaic mammals
summary(lm3)  #*********

lm4<-glm(Faunal_Feve~+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm4)
plot(lm4)

lm5<-glm(Faunal_Feve~MixedModel_Data$lognormal+MixedModel_Data$do18+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm5<-glm(Faunal_Feve~MixedModel_Data$lognormal+MixedModel_Data$do18, data=MixedModel_Data)
summary(lm5)
plot(lm5)

lm1<-glm(Faunal_Feve~MixedModel_Data$lognormal+MixedModel_Data$Faunal_Prop.Archaic, data=MixedModel_Data) #species richness lognormal
summary(lm1)




##################################################################################################3
#Null Model trying to pull 10 random communities from  Paleocene and then from the first 10MA of Eocene
Communities<-read.csv("FunctionalCommunities_Final.csv")
Eocene<-subset(Communities, Communities$Ma <= 56 & Communities$Ma >= 46)
Paleocene<-subset(Communities, Communities$Ma >= 56)

rownames(Eocene)<-(1:5)#Setting rownames to 1-29 to link with numbers randomized below

Null <- replicate(100,sample(unique(Eocene$Feve), 5))  #This takes the Fric directly so I don't have to match
Null<-write.csv(Null, file='Paleocene_NullModelFric.csv')
Null1<-read.csv('Feve_NullModel.csv')

#FRIC
#Plotting in Histogram
ggplot(Null1, aes(Eocene_Mean)) +
  geom_histogram(bins=15,fill='khaki1', color='black')+
 xlim(0.585, 0.91)+
  geom_vline(aes(xintercept=mean(Eocene_Mean)),
             color="blue", linetype="dashed", size=3)+
 # geom_vline(aes(xintercept=mean(Null1$Eocene_Mean)),
  #           color="orange1", linetype="dashed", size=2)+
  theme_classic()


t.test(Null1$Eocene_Mean, Null1$Paleocene_Mean)



###############################################################################################
#Run a breakpoint analysis and include loess regression
Comm<-read.csv("FunctionalCommunities_Final.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Fric", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFric<-BetaD[!duplicated(BetaD), ] 
na.omit(BetaFric)


Comm<-read.csv("FunctionalCommunities_Final.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Fdiv", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFdiv<-BetaD[!duplicated(BetaD), ] 
na.omit(BetaFdiv)



Comm<-read.csv("FunctionalCommunities_Final.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Feve", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFeve<-BetaD[!duplicated(BetaD), ] 
na.omit(BetaFeve)


#Find best span for each metric
library(fANCOVA)
FTSE.Fric <- loess.as(BetaFric$Fric, BetaFric$Ma, degree = 2, criterion = c("aicc","gcv")[1], user.span = 0.35, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Fric, data.frame(Index=Comm$Ma))

FTSE.Fdiv <- loess.as(BetaFdiv$Fdiv,BetaFdiv$Ma, degree = 2, criterion = c("aicc","gcv")[1], user.span = 0.35, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Fdiv, data.frame(Index=Comm$Ma))

FTSE.Feve <- loess.as(BetaFeve$Feve, BetaFeve$Ma, degree = 2, criterion = c("aicc","gcv")[1], user.span = 0.35, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Feve, data.frame(Index=Comm$Ma))


#####################################################################################
#Run a breakpoint analysis and include loess regression
Comm<-read.csv("1myr_TimeBin_Data_Mixed_Model.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Faunal_Fric", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFFric<-BetaD[!duplicated(BetaD), ] 
BetaFFric<-na.omit(BetaFFric)


Comm<-read.csv("1myr_TimeBin_Data_Mixed_Model.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Faunal_Fdiv", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFFdiv<-BetaD[!duplicated(BetaD), ] 
BetaFFdiv<-na.omit(BetaFFdiv)



Comm<-read.csv("1myr_TimeBin_Data_Mixed_Model.csv")
na.omit(Comm)
Comm<-Comm[, (names(Comm) %in% c("Faunal_Feve", "Ma"))] #do for each FD indice
Comm<-as.data.frame(Comm)
BetaD<- Comm[order(Comm$Ma),]
BetaFFeve<-BetaD[!duplicated(BetaD), ] 
BetaFFeve<-na.omit(BetaFFeve)


#Find best span for each metric
library(fANCOVA)
FTSE.Fric <- loess.as(BetaFFric$Faunal_Fric, BetaFFric$Ma, degree = 1, criterion = c("aicc","gcv")[1], user.span = 0.5, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Fric, data.frame(Index=Comm$Ma))

FTSE.Fdiv <- loess.as(BetaFFdiv$Faunal_Fdiv,BetaFFdiv$Ma, degree = 1, criterion = c("aicc","gcv")[1], user.span = 0.5, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Fdiv, data.frame(Index=Comm$Ma))

FTSE.Feve <- loess.as(BetaFFeve$Faunal_Feve, BetaFFeve$Ma, degree = 1, criterion = c("aicc","gcv")[1], user.span = 0.5, plot = T)
FTSE.lo.predict3 <- predict(FTSE.Feve, data.frame(Index=Comm$Ma))







#########################################################################
#Breakpoint Analysis
#Run for each metric from each timescale

#Determine the number of appropriate breakpoints
BetaFB<-BetaD[, (names(BetaFric) %in% ("Fric"))]
BetaFB<-as.data.frame(BetaFB)

library(breakpoint)
CE.Normal.MeanVar(BetaFB, Nmax = 4, eps = 0.05, rho = 0.05, M = 250, h = 5, a = 0.8, b = 0.8,
                  distyp = 1, penalty = "AIC", parallel = FALSE)
#Breakpoints analysis
#Using segemented to create a linear model and determine breakpoints
# create a linear model
my.lm <- lm(Fric~Ma, data = BetaFric) 
summary(my.lm)
# Extract te coefficients from the overall model
my.coef <- coef(my.lm)

# add the regression line to the graph
# setting the aesthetics to a constant - this provides a name that we can reference later when we add additional layers
p<- geom_abline(intercept = my.coef[1], 
                slope = my.coef[2])
plot(my.lm)

# have to provide estimates for breakpoints from above
library(segmented)
my.seg <- segmented(my.lm, seg.Z = ~ Ma, psi = list(Ma =c(58)))

# display the summary
summary(my.seg)
library(wiqid)
AICc(my.seg)

Fbsorsummary<-my.seg$psi
slope(my.seg)

# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(Distance = BetaFric$Bin, Elevation = BetaFric$Fric)
my.Fbsor<- my.seg$psi[, 2]
my.lines <- my.seg$psi[, 2]

# plot the fitted model
ggplot(my.model, aes(x = Distance, y = Elevation)) + 
  geom_smooth()+
  geom_vline(xintercept = my.lines, linetype = "solid")


