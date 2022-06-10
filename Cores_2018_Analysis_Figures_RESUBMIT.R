# Project: Variation in Sediment Phosphorus (P) Manuscript - JGR Biogeosciences RESUBMISSION (analysis updates)
# Last Modified 10 June 2022
# Contributers: Ellen Albright, Dr. Grace Wilkinson
# Description: The following script provides all analyses and figures for our manuscript "High inter- and intra-lake variation in sediment phosphorus pools in shallow lakes"

# Data citation: Albright, Ellen; Wilkinson, Grace; Fleck, Rachel; Shingai, Quin (2020): Analysis of Sediment Phosphorus Pools in Shallow Lakes. figshare. Dataset. https://doi.org/10.6084/m9.figshare.13362971.v1 

# The code can be run with the following datasets (NOTE: for each dataset, there is a README file with detailed meta-data)
  # "Cores_ALLDATA_2018.csv" - sediment P chemistry and physical characteristics for 7 shallow lakes in NW Iowa, USA. 
  # "Chla_MobileP_Regression.csv" - summarizes mobile (redox-sensitive, loosely-bound, organic) sediment P fractions and long term chlorophyll a concentrations for each lake
  # "Chla_MobileP_Regression_NEW.csv" - summarizes mobile P in 0-6cm slice, and long-term Chla through 2018
  # "PCA_envfit.csv" - contains lake and watershed variables as potential interpretations of PCA biplot of profundal sediment P speciation
  # "Swan_macrophytes_2018.csv" - macrophyte survey across Swan Lake (need "studylakes2018.shp" if you want to map the lake outline)
  # "Swan_SedP_2018.csv" - summary of sediment P results and sampling site coordinates for mapping
  # "Rarefaction_RMSE.csv" - summary of RMSE values for TP and loosely-bound P concentrations from a rarefaction analysis based on number of sampling sites
  # "Dv_Analysis.csv" - includes morphometric data to calculate lake volume development. Used for Supporting Information Figure.

# Clear environment, set working directory - UPDATE AS NEEDED ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
getwd()
setwd('C:/Users/Ellen/Desktop/Box Sync/Albright DISSERTATION/Albright_Chapter_1/Data/Data_Scripts_SUBMIT')

# REQUIRED PACKAGES FOR CODE - install as needed
# install.packages("tidyverse")
# install.packages("RColorBrewer")
# install.packages("robCompositions")
# install.packages("vegan")
# install.packages('lmerTest')
# install.packages("sf")
# install.pacakges("ggnewscale")
# install.packages("ggnewscale")
# install.packages("infer")
library(tidyverse) #general, data cleaning
library(RColorBrewer) #general, data visualization
library(robCompositions) #analysis, CoDA
library(vegan) #analysis, PCA
library(lmerTest) #analysis, mixed model
library(sf) #analysis, spatial
library(ggnewscale) #visualization, spatial
library(gridExtra) #general, visualization
library(infer) #rarefaction analysis, resampling

# COLORBREWER PALLETES, etc.
par(mfrow=c(1,1))
display.brewer.all()

### OVERVIEW OF SCRIPT ORGANIZATION ### -----------------------------------------------------------------------------------------------------------------------------------------------
# PART 1: Descriptions of inter-lake variation in profundal sediment P
#          1a.calculations of various summary statistics (Table 2) 
# PART 2: Analysis of inter-lake variation in profundal sediment P
#          2a. compositional data analysis of profundal sediment P speciation among lakes (Figure 1A-B)
# PART 3: Further analysis of water quality implications of inter-lake variation in profundal sediment P
#          3a. simple linear regression of long-term chlorophyll a concentrations as a function of mobile sediment P species (Figure 2)
# PART 4: Intra-lake variation in sediment P - horizontal variation across the lakebed 
#          4a. mixed model regression effects of water depth at the coring location on total sediment P and loosely-bound sediment P by lake (Figure 3)
#          4b. calculations of coefficient of variation for total and loosely-bound P (Table 2)
# PART 5: Intra-lake variation in sediment P - Macrophyte beds and sediment P Pools
#          5a. Within-lake variation in sediment loosely-bound and total P in Swan Lake (Figure S2)
# PART 6: Intra-lake variation in sediment P - Rarefaction analysis to determine minimum sampling frequency needed
#          6a. Rarefaction curves for total and loosely-bound P for each study lake (Figure 5A-B)
# PART 7: Intra-lake variation in sediment P and basin volume development (Figure 4)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### PART 1 - DESCRIPTIONS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (TABLE 2) ----------------------------------------------------------------------------------------------------
# Read in full 2018 sediment dataset and subset out samples from the deep site of the study lakes (SiteType=="Deep")
down<-read.csv("Cores_ALLDATA_2018.csv")
deep<-subset(down,SiteType=="Deep")

### PART 1a. Calculations of summary statistics for Table 2 ----------------------------------------------------------------------------------------------------------------------------
# Mean and standard deviation of total P over all intervals in the deep site core for each lake
tapply(deep$TP_ug,deep$LakeName,mean,na.rm=TRUE)
tapply(deep$TP_ug,deep$LakeName,sd,na.rm=TRUE)
tapply(deep$LooseP,deep$LakeName,mean,na.rm=TRUE)

# Average percent contribution of each sediment P fraction
# Need to tidy the dataset a bit before making these calculations because we do not have labile organic P data for all of our lakes
result<-deep %>% 
  filter(SampleID != "CE1B" & SampleID != "CE1F" & SampleID != "NT1A") %>%  #filter out rows that have missing TP data
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,LabileOP,CalciumP,TP_ug) #select relevant columns to work with
result1<-result %>%
  filter(LakeID %in% c("FI","CE","SO","SI")) %>% #first working with 4 lakes that have labile organic P data
  mutate(RefractoryP = TP_ug - (LooseP+IronP+AlP+CalciumP+LabileOP), #calculate refractory organic P concentration and then relative proportions of each fraction
         Per_Loose=LooseP/TP_ug,Per_Fe=IronP/TP_ug,Per_Al=AlP/TP_ug,Per_Ca=CalciumP/TP_ug,
         Per_Labile=LabileOP/TP_ug,Per_Refrac=RefractoryP/TP_ug,Per_Organic=(RefractoryP+LabileOP)/TP_ug,
         Per_Mobile=Per_Loose+Per_Fe+Per_Organic, Per_IronLoose =Per_Loose+Per_Fe )%>%
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,TP_ug,
         Per_Loose,Per_Fe,Per_Al,Per_Ca,Per_Labile,Per_Refrac,Per_Organic,Per_Mobile,Per_IronLoose) #add refractory P if you want to calculate percent of total!

tapply(result1$Per_Loose,result1$LakeName,mean) # percent loosely-bound P, averaged across all intervals in the deep site core for each lake
tapply(result1$Per_Fe,result1$LakeName,mean) # percent redox-sensitive P
tapply(result1$Per_Al,result1$LakeName,mean) # percent aluminum-bound P
tapply(result1$Per_Ca,result1$LakeName,mean) # percent calcium-bound P
tapply(result1$Per_Labile,result1$LakeName,mean) # percent labile organic P
tapply(result1$Per_Refrac,result1$LakeName,mean) # percent refractory organic P
tapply(result1$Per_Organic,result1$LakeName,mean) # percent total organic P (as sum of labile and refractory)

# NUBMERS USED FOR CHLA REGRESSION - MEAN AND SD of MOBILE P OVER DEEP CORE
tapply(result1$Per_Mobile,result1$LakeName,mean) # percent mobile P, averaged across all intervals in the deep site core for each lake
tapply(result1$Per_Mobile,result1$LakeName,sd) # percent mobile P, sd across all intervals in the deep site core for each lake

# Average percent contribution of each sediment P fraction - lakes without labile P data
result2<-result %>%
  filter(LakeID %in% c("NT","ST","SW")) %>% #next working with 3 lakes that do NOT have labile organic P data
  mutate(OrganicP = TP_ug - (LooseP+IronP+AlP+CalciumP), # total organic P concentration and then relative proportions of each fraction
         Per_Loose=LooseP/TP_ug,Per_Fe=IronP/TP_ug,Per_Al=AlP/TP_ug,Per_Ca=CalciumP/TP_ug,Per_Organic=OrganicP/TP_ug,
         Per_Mobile=Per_Loose+Per_Fe+Per_Organic, Per_IronLoose =Per_Loose+Per_Fe) %>%
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,TP_ug,
         Per_Loose,Per_Fe,Per_Al,Per_Ca,Per_Organic,Per_Mobile,Per_IronLoose) #add refractory P if you want to calculate percent of total!

tapply(result2$Per_Loose,result2$LakeName,mean) # percent loosely-bound P, averaged across all intervals in the deep site core for each lake
tapply(result2$Per_Fe,result2$LakeName,mean) # percent redox-sensitive P
tapply(result2$Per_Al,result2$LakeName,mean) # percent aluminum-bound P
tapply(result2$Per_Ca,result2$LakeName,mean) # percent calcium-bound P
tapply(result2$Per_Organic,result2$LakeName,mean) # percent total organic P (assumed to be difference between total P and 4 measured fractions)

# NUBMERS USED FOR CHLA REGRESSION - MEAN AND SD of MOBILE P OVER DEEP CORE
tapply(result2$Per_Mobile,result2$LakeName,mean) # percent mobile P, averaged across all intervals in the deep site core for each lake
tapply(result2$Per_Mobile,result2$LakeName,sd) # percent mobile P, sd across all intervals in the deep site core for each lake

### PART 2 - ANALYSIS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (TABLE 2, FIGURE 1A-B) -------------------------------------------------------------------------------------------
### PART 2a. compositional data analysis (CoDA) of profundal sediment P speciation among lakes (Figure 1A-B) -----------------------------------------------------------------------------
# Uses "deep" dataframe from Part 1

# First format data for compositional data analysis
deep2<-deep %>% 
  filter(SampleID != "CE1B" & SampleID != "CE1F" & SampleID != "NT1A") %>%  #filter out rows that have missing TP data 
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,LabileOP,CalciumP,TP_ug) #select relevant columns to work with

# NOTE: the CoDA is based on the concentrations of loosely-bound, redox-sensitive, calcium-bound, aluminum-bound, and total organic P
#       total organic P is the sum of labile and refractory components; therefore, we needed to work with lakes that didn't have labile organic P measurements separately

deep3<-deep2 %>%
  filter(LakeID %in% c("FI","CE","SO","SI")) %>% #first working with 4 lakes that have labile organic P data
  mutate(RefractoryP = TP_ug - (LooseP+IronP+AlP+CalciumP+LabileOP), #calculate refractory organic P concentration and then relative proportions of each fraction
         OrganicP=RefractoryP+LabileOP) %>% #combining labile and refractory organic P proportions to have a variable that applies to other three lakes
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,OrganicP,TP_ug) 
deep4<-deep2 %>%
  filter(LakeID %in% c("SW","NT","ST")) %>% #3 lakes without labile OP data
  mutate(OrganicP = TP_ug - (LooseP+IronP+AlP+CalciumP)) %>% #calculate total organic P concentration and then relative proportions of each fraction
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,OrganicP,TP_ug)
core<-full_join(deep4,deep3)
View(core) #the concentrations of the 5 P fractions are in columns 6-10

# Time for the CoDA - need the RobCompositions package
# Take the centered log raio (CLR) of the compositional data
CLR<-cenLR(core[, 6:10])

# the cenLR function returns two lists: x.clr gives the CLR of each concentration, gm gives the geometric means for each row (sample)
# for now, we want the CLR values in a data frame so that we can add these values to the main dataset and then perform a PCA
CLRdf<-data.frame(CLR$x.clr)

# Join CLR dataframe to "core" dataframe - need to change CLRdf column names first
CLRdf2<-rename(CLRdf,CLR_loose=LooseP,CLR_redox=IronP,CLR_Al=AlP,CLR_Ca=CalciumP,CLR_Organic=OrganicP)
core_CLR<-cbind(core,CLRdf2)
View(core_CLR)# CLR values in columns 12-16

# PCA on covariance matrix using rda() function, need vegan package
core.rda<-rda(core_CLR[,12:16],scale=FALSE) #selects all rows but only columns 12-16, scale=FALSE uses covariance matrix
summary(core.rda) # PC1 explains 45.27% of variation, PC2 explains 36.42% (sum=81.69%)

# PCA Bioplot (Figure 1A-B)
with(core,levels(LakeName))
brewer.pal(7,"Dark2")
CEcol<-"#A6761D" #brown
FIcol<-"#D95F02" #orange
NTcol<-"#7570B3" #purpleblue
SIcol<-"#E7298A" #pink
SOcol<-"#E6AB02" #mustard
STcol<- "#1B9E77" #teal
SWcol<-"#66A61E" #green
colvec<-c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol)

### Adding some environmental data to PCA biplot
env <- read.csv("PCA_envfit.csv")
View(env)

# envfit() function to add environmental data 
x_clr<-envfit(core.rda,core_CLR[,12:16],choices=c(1,2)) #"species" scores
x_landcover<-envfit(core.rda,env[,c(17:18,20:21)],choices=c(1,2)) #landcover
x_env<-envfit(core.rda,env[,c(10,12,15,16)],choices=c(1,2)) #lake variables

# Set up Figure 1A-B
png("Figure1.png", width=5, height=9, units="in",res=300)
par(mfrow=c(2,1), mai=c(0.8,0.8,0.1,0.1), omi=c(0.4,0.05,0,0))

# 1A: PCA biplot with "species" scores (CLR-transformed sediment P concentrations)
biplot(core.rda,choices=c(1,2),display="sites",type="points",scaling=3,xlab="",ylab="",ylim=c(-1,1),xlim=c(-1,1)) #plot PCA, symmetrical scaling, "sites" only
title(xlab="PC1 (45.27)%",ylab="PC2 (36.42%)",cex.lab=1,) #show % variation explained by each principal component
ordihull(core.rda,group=core$LakeName,scaling=3,label=FALSE,draw="polygon",col=c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol),border="white",alpha=0.3)
with(core,points(core.rda,display="sites",col="white",pch=21,scaling=3,bg=colvec[LakeName],cex=1.4)) #color points by lake
text(-0.65,-0.30,labels="Swan",col="#3e6512",cex=1)
text(0.29,-0.69,labels="North Twin",col="#4f4a8c",cex=1)
text(0.68, -0.2,labels="South Twin",col="#105d46",cex=1)
text(0.55,0.1,labels="Silver",col="#b01463",cex=1)
text(0.16,0.6,labels="Storm",col="#9a7301",cex=1)
text(-0.25,0.28,labels="Center",col="#654812",cex=1)
text(-0.31,-0.05,labels="Five Island",col="#8d3e01",cex=1)
plot(x_clr,choices=c(1,2),at=c(0,0),axis=FALSE,p.max=NULL,add=TRUE,col="black",labels=c(" "," "," "," "," "),cex=1.5)
text(-0.72,0.44,labels="Redox-Sensitive P",cex=1)
text(-0.99,0.18,labels="Al-Bound
     P",cex=1)
text(0.05,-0.88,labels="Organic P",cex=1)
text(0.7,-0.52,labels="Loosely-Bound P",cex=1)
text(0.83,0.72,labels="Ca-Bound P",cex=1)
text(-1,0.95,labels="A",cex=1.5)

# Figure 1B - All environmental data
biplot(core.rda,choices=c(1,2),display="sites",type="points",scaling=3,xlab="",ylab="",ylim=c(-1,1),xlim=c(-1,1)) #plot PCA, symmetrical scaling, "sites" only
title(xlab="PC1 (45.27)%",ylab="PC2 (36.42%)",cex.lab=1,) #show % variation explained by each principal component
ordihull(core.rda,group=core$LakeName,scaling=3,label=FALSE,draw="polygon",col="gray50",border="white",alpha=0.3)
with(core,points(core.rda,display="sites",col="white",pch=21,scaling=3,bg="gray60",cex=1.5)) 
plot(x_landcover,choices=c(1,2),at=c(0,0),axis=FALSE,p.max=NULL,add=TRUE,col="darkorange2",labels=c(" "," "," "," "),cex=1.5) #land cover!
text(-0.8,0.45,labels="Grassland",cex=1,col="darkorange2")
text(-0.95,0.22,labels="Forest",cex=1,col="darkorange2")
text(0.55,0.2,labels="Developed",cex=1,col="darkorange2")
text(0.8,-0.5,labels="Cropland",cex=1,col="darkorange2")
plot(x_env,choices=c(1,2),at=c(0,0),axis=FALSE,p.max=NULL,add=TRUE,col="dodgerblue2",labels=c(" "," "," "," "),cex=1.5) #lake variables
text(-0.44,0.7,labels="Maximum Depth",cex=1,col="dodgerblue2")
text(0.4,0.84,labels="Bulk Density",cex=1,col="dodgerblue2")
text(0.7,-0.64,labels="Volume Development",cex=1,col="dodgerblue2")
text(0.6,-0.89,labels="Sediment Organic Matter",cex=1,col="dodgerblue2")
text(-1,0.95,labels="B",cex=1.5)
dev.off()

# From Audrey: "I ran the envfit function on both the untransformed and the clr-transformed data, and the plots were really similar.
#   It would be better to use the clr-transformed data since that's what you plugged into the rda function. 
#   But because they're really close, I think you could report the untransformed data and say that you ran it on the clr-transformed data and the results didn't change.


### PART 3 - FURTHER ANALYSIS OF WATER QUALITY IMPLICATIONS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (FIGURE 2) ----------------------------------------------------------------
### PART 3a. simple linear regression of long-term chlorophyll a concentrations as a function of mobile sediment P species (Figure 2) -------------------------------------------------

###      Step 1: re-clean the ALM Chlorophyll data to exclude 2019
chla<-read.csv("chlorophyll_AQUIA.csv")

chla<-chla %>% 
  select(name, sampleDate, year,result,unit) %>% 
  filter(year!="2019")

# calculate mean, sd, and n for each lake
tapply(chla$result, chla$name, mean)
tapply(chla$result, chla$name, sd)
length(unique(chla$result[chla$name=="Center Lake"]))
length(unique(chla$result[chla$name=="Five Island Lake"]))
length(unique(chla$result[chla$name=="North Twin Lake (maximum water depth)"]))
length(unique(chla$result[chla$name=="Silver Lake Max Depth"]))
length(unique(chla$result[chla$name=="South Twin Lake"]))
length(unique(chla$result[chla$name=="Storm Lake Max Depth"]))
length(unique(chla$result[chla$name=="Swan Lake Max Depth"]))

###      Step 2: re-calculate mobile P fractions for the 0-6cm slice only!
core<-read.csv("Cores_ALLDATA_2018.csv")
core<-core %>% 
  filter(SiteType=="Deep") %>% 
  select(LakeName,SampleID,TopDepth,LooseP,IronP,AlP,LabileOP,CalciumP,TP_ug) %>% 
  drop_na(TP_ug) %>% 
  filter(TopDepth<=4)

core$mobile<-core$TP_ug-(core$CalciumP+core$AlP)
core$mobilePER<-(core$mobile/core$TP_ug)*100

tapply(core$mobilePER, core$LakeName, mean)
tapply(core$mobilePER, core$LakeName, sd)

tapply(core$mobile, core$LakeName, mean)
tapply(core$mobile, core$LakeName, sd)


###      Step 3: Read in auxiliary dataset that summarizes mobile sediment P species and long-term chla concentrations
chla<-read.csv("Chla_MobileP_Regression_NEW.csv")
View(chla)

# Define variables
x<-chla$mobilePercent4
SDx<-(chla$mobilePERSD4/sqrt(chla$mobileN4))
xlabel<-"% Mobile P of Total Sediment P"

y<-chla$chlorophyll
SDy<-chla$chlorophyllSD/sqrt(chla$chlorophyllN)

# simple linear model of chla as a function of mobile sediment P contribution
mod<-lm(y~x)
summary(mod) #Adj R2=0.5232, F(1,5)=7.584, p=0.04013
#confidence interval of the estimated B1 parameter (slope):
confint(mod,level=0.95) # Fitted B1 is 3.953 with a 95% confidence interval of (0.2630499, 7.642428)

Chlamodel<-ggplot(data=chla,aes(x=x,y=y)) +
  geom_smooth(method="lm",formula=y~x,se=TRUE,level=0.95,color='black',alpha=0.1) +
  geom_point(size=2,color="cadetblue4",alpha=0.8) +  
  geom_errorbar(aes(ymin=y-SDy, ymax=y+SDy),width=0,color="cadetblue4")+geom_errorbar(aes(xmin=x-SDx, xmax=x+SDx),width=0,color="cadetblue4")+
  geom_point(size=2,color="cadetblue4")+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlim(80,100) + ylim(0,110) + xlab("Mobile P (% Total Sediment P)") + 
  ylab(bquote('Long-Term Mean Chlorophyll (mg '  *L^-1*')'))+
  theme(axis.text=element_text(color="black",size=10),axis.title=element_text(size=11))
ggsave("Figure2.png",Chlamodel,width=3,height=3,units="in",dpi=300)

### PART 4 - INTRA-LAKE VARIATION IN SEDIMENT P - HORIZONTAL VARIATION ACROSS THE LAKEBED (FIGURE 3, TABLE 2) -------------------------------------------------------------------------
### PART 4a. Mixed model regression effects of water depth at the coring location on total sediment P and loosely-bound sediment P by lake (Figure 3) ---------------------------------

# This analysis uses the "down" dataframe originally read in for Part 1A.
# We will need to format this dataframe to select all of the "littoral" site cores and the average value of all the sediment slices of the deep site core for each lake
# First - subset the whole dataframe ("down") to isolate the littoral sites, then select only the columns we need
lit<-subset(down,SiteType=="Littoral")
lit2<-lit %>%
  select(LakeName,WaterDepth,LooseP,TP_ug)

# Next - take the average values over the core interval slices for the deep site cores for each lake
deep_av1<-deep %>% 
  filter(SampleID != "CE1B" & SampleID != "CE1F" & SampleID != "NT1A")  #filter out rows that have missing TP data
deep_av<-aggregate(cbind(WaterDepth,LooseP,TP_ug)~LakeName,deep_av1,mean,na.rm=TRUE)

# Finally - bind the littoral and averaged deep site dataframe together
mix<-rbind(lit2,deep_av)

# Time for the mixed models! Need lmerTest package loaded
# Loosely-bound P Mixed Effects model====================================================
looseP_mixed<-lmerTest::lmer(LooseP ~ WaterDepth + (1|LakeName), data=mix, na.action="na.fail", REML=F)
summary(looseP_mixed)
confint(looseP_mixed,level=0.95) #B1 water depth = 12.787 [7.180, 18.334]
#Calculate the variance explained by the model (R-squared)
rsq<-function(emp,pred){
  sstot<-sum((emp-mean(emp))^2)
  ssres<-sum((pred-emp)^2)
  return(1-(ssres/sstot))
}
rsq(mix$LooseP, predict(looseP_mixed)) #R2=0.846
# Null model - don't allow intercept to vary by lake
looseP_null<-lm(LooseP ~ WaterDepth,data=mix)
summary(looseP_null)
#Compare the mixed model to the null model
anova(looseP_mixed, looseP_null) #including random effects on model intercept improved fit (p<0.0001)

#Total P Mixed Effects Model====================================================
totalP_mixed<-lmerTest::lmer(TP_ug ~ WaterDepth + (1|LakeName), data=mix, REML=F)
summary(totalP_mixed)
confint(totalP_mixed,level=0.95) #B1 water depth = 65.28 [27.5661, 102.4886]

#Calculate the variance explained by the model (R-squared)
rsq(mix$TP_ug, predict(totalP_mixed)) #R2=0.631
# Null model - don't allow lake to be a random effect on the intercept
totalP_null<-lm(TP_ug ~ WaterDepth,data=mix)
summary(totalP_null)
#Compare the mixed model to the null model
anova(totalP_mixed, totalP_null) #including random effects on model intercept improved fit (p<0.0001)

# PLOT THE DATA (Figure 3) ============================================================
png("Figure3.png", width=4, height=6.5, units="in",res=300)
par(mfrow=c(2,1), mai=c(0.5,0.75,0.1,0.1), omi=c(0.4,0.05,0,0))

#Build the color palette for plotting
pttrans=0.5 #transparency value for points
ptpal=brewer.pal(7,"Dark2")
ptrgb<-col2rgb(ptpal)
ptpal<-rgb(red=ptrgb[1,],green=ptrgb[2,],blue=ptrgb[3,],alpha=pttrans*255, maxColorValue = 255) #color palette for transparent points
lnpal=brewer.pal(7,"Dark2") #color palette for solid lines

plot(NA, NA, xlim=c(0,7), ylim=c(0,1400), 
     xlab="", ylab="", cex.lab=1, cex.axis=0.85)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey99",border=NA)
mtext(side=2, line=2.5, expression(Total~P~"("*mu*g~P~g^-1~dry~sed*")"))
#Center
points(mix[mix$LakeName=="Center", "WaterDepth"], mix[mix$LakeName=="Center", "TP_ug"],pch=16, col=ptpal[7], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[1,1],b=coef(totalP_mixed)$LakeName[1,2], col=lnpal[7], lwd=2)
#Five Island
points(mix[mix$LakeName=="Five Island", "WaterDepth"], mix[mix$LakeName=="Five Island", "TP_ug"],pch=16, col=ptpal[2], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[2,1],b=coef(totalP_mixed)$LakeName[2,2], col=lnpal[2], lwd=2)
#North TWin
points(mix[mix$LakeName=="North Twin", "WaterDepth"], mix[mix$LakeName=="North Twin", "TP_ug"],pch=16, col=ptpal[3], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[3,1],b=coef(totalP_mixed)$LakeName[3,2], col=lnpal[3], lwd=2)
#Silver
points(mix[mix$LakeName=="Silver", "WaterDepth"], mix[mix$LakeName=="Silver", "TP_ug"],pch=16, col=ptpal[4], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[4,1],b=coef(totalP_mixed)$LakeName[4,2], col=lnpal[4], lwd=2)
#South Twin
points(mix[mix$LakeName=="South Twin", "WaterDepth"], mix[mix$LakeName=="South Twin", "TP_ug"],pch=16, col=ptpal[1], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[5,1],b=coef(totalP_mixed)$LakeName[5,2], col=lnpal[1], lwd=2)
#Storm
points(mix[mix$LakeName=="Storm", "WaterDepth"], mix[mix$LakeName=="Storm", "TP_ug"],pch=16, col=ptpal[6], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[6,1],b=coef(totalP_mixed)$LakeName[6,2], col=lnpal[6], lwd=2)
#Swan
points(mix[mix$LakeName=="Swan", "WaterDepth"], mix[mix$LakeName=="Swan", "TP_ug"],pch=16, col=ptpal[5], cex=1.25)
abline(a=coef(totalP_mixed)$LakeName[7,1],b=coef(totalP_mixed)$LakeName[7,2], col=lnpal[5], lwd=2)
#LEGEND
legend("bottomright", legend=c("Center", "Five Island", "North Twin", "Silver", "South Twin", "Storm", "Swan"),
       col = c(lnpal[7], lnpal[2], lnpal[3], lnpal[4], lnpal[1], lnpal[6], lnpal[5]), pch=15, pt.cex=1.5, cex=0.8, ncol=2)
box()

plot(NA, NA, xlim=c(0,7), ylim=c(0,250), 
     xlab="", ylab="", cex.lab=1, cex.axis=0.85)
rect(xleft=-1,ybottom=-1,xright=60,ytop=0,col="grey99",border=NA)
mtext(side=1, line=2.5, "Water Depth at Coring Location (m)")
mtext(side=2, line=2.5, expression(Loosely-bound~P~"("*mu*g~P~g^-1~dry~sed*")"))
#Center
points(mix[mix$LakeName=="Center", "WaterDepth"], mix[mix$LakeName=="Center", "LooseP"],pch=16, col=ptpal[7], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[1,1],b=coef(looseP_mixed)$LakeName[1,2], col=lnpal[7], lwd=2)
#Five Island
points(mix[mix$LakeName=="Five Island", "WaterDepth"], mix[mixh$LakeName=="Five Island", "LooseP"],pch=16, col=ptpal[2], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[2,1],b=coef(looseP_mixed)$LakeName[2,2], col=lnpal[2], lwd=2)
#North TWin
points(mix[mix$LakeName=="North Twin", "WaterDepth"], mix[mix$LakeName=="North Twin", "LooseP"],pch=16, col=ptpal[3], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[3,1],b=coef(looseP_mixed)$LakeName[3,2], col=lnpal[3], lwd=2)
#Silver
points(mix[mix$LakeName=="Silver", "WaterDepth"], mix[mix$LakeName=="Silver", "LooseP"],pch=16, col=ptpal[4], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[4,1],b=coef(looseP_mixed)$LakeName[4,2], col=lnpal[4], lwd=2)
#South Twin
points(mix[mix$LakeName=="South Twin", "WaterDepth"], mix[mix$LakeName=="South Twin", "LooseP"],pch=16, col=ptpal[1], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[5,1],b=coef(looseP_mixed)$LakeName[5,2], col=lnpal[1], lwd=2)
#Storm
points(mix[mix$LakeName=="Storm", "WaterDepth"], mix[mix$LakeName=="Storm", "LooseP"],pch=16, col=ptpal[6], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[6,1],b=coef(looseP_mixed)$LakeName[6,2], col=lnpal[6], lwd=2)
#Swan
points(mix[mix$LakeName=="Swan", "WaterDepth"], mix[mix$LakeName=="Swan", "LooseP"],pch=16, col=ptpal[5], cex=1.25)
abline(a=coef(looseP_mixed)$LakeName[7,1],b=coef(looseP_mixed)$LakeName[7,2], col=lnpal[5], lwd=2)

dev.off()

### PART 4b. Calculations of coefficient of variation for total and loosely-bound P (Table 2) -----------------------------------------------------------------------------------------
meanTP<-tapply(mix$TP_ug,mix$LakeName,mean)
sdTP<-tapply(mix$TP_ug,mix$LakeName,sd)
CVTP<-(sdTP/meanTP)*100
CVTP

meanloose<-tapply(mix$LooseP,mix$LakeName,mean)
sdloose<-tapply(mix$LooseP,mix$LakeName,sd)
CVloose<-(sdloose/meanloose)*100
CVloose

### PART 5 - INTRA-LAKE VARIATION IN SEDIMENT P - MACROPHYTE BEDS AND SEDIMENT P POOLS  (Figure S2) --------------------------------------------------------------------
### PART 5a. Within-lake variation in sediment loosely-bound and total P in Swan Lake (Figure S2)---------------------------------------------------------------------------------------

# Read in shapefiles for lake outline (make sure sf package is loaded)
lakes<-st_read("studylakes2018.shp") #old shape file from MC project- has outlines of all ALM lakes
st_crs(lakes) #meta-data on shape file - coordinate system and projection
names(lakes) #column names for shp attribute table
levels(lakes$Lake_Name) #tells you all the values in the lake name colume (get the correct spelling to select specific lake)
swan_shp<-lakes %>% #select the swan lake outline only
  filter(Lake_Name=="Swan Lake")

# Read in SW macrophyte survey data
SWmac<-read.csv("Swan_macrophytes_2018.csv")
SWmac<-SWmac %>%
  filter(TOTAL!="NA")
# Read in SW sediment core data
SWcore<-read.csv("Swan_SedP_2018.csv")

# Convert both to spatial files
sites<-st_as_sf(SWmac,coords=c("LONGUSE","LATUSE"),crs=4893,agr="constant") #convert dataframe with Swan macrophyte sites to a spatial file (4893 is the CRS code for NAD83 - what swan outline shp projected in)
cores<-st_as_sf(SWcore,coords=c("LONG","LAT"),crs=4893,agr="constant")
names(sites)
names(cores)

# Plot data together (Figure S1)
SW_plot<-ggplot() + 
  geom_sf(data = swan_shp, size = 1, color = "black", fill = "white") + #f2f8f9 #e0eff1
  geom_sf(data=sites,aes(size=TOTAL),col="#b8e186")+scale_size(range=c(5,11))+
  new_scale("size")+
  geom_sf(data=cores,aes(size=TP_ug,fill=LoosePER),shape=21,col="black")+scale_size(range=c(5,16))+scale_fill_continuous(low=("#fde0ef"),high=("#8e0152"))+
  theme_void()+theme(legend.position="none")+
  coord_sf()
ggsave("Figure4.png",SW_plot,width=5,height=4,units="in",dpi=300)


### PART 6 - Implications of horizontal variation in sediment P pools --------------------------------------------------------------------------------------------
### Rarefaction analysis to determine minimal sampling needed
#        6a. Rarefaction curves for total and loosely-bound P for each study lake (Figure 5A-B)

# First, need to wrangle the data a bit - re-doing the data cleaning in PART 4a. 
# Goal - df of all the littoral sites AND the average over the slices from the deep core (so 10 sites per lake)
# First - subset the whole dataframe ("down") to isolate the littoral vs. deep sites, then select only the columns we need
deep<-subset(down,SiteType=="Deep")
lit<-subset(down,SiteType=="Littoral")
lit2<-lit %>%
  select(LakeName,LooseP,TP_ug)

# Next - take the average values over the core interval slices for the deep site cores for each lake
deep_av1<-deep %>% 
  filter(SampleID != "CE1B" & SampleID != "CE1F" & SampleID != "NT1A")  #filter out rows that have missing TP data
deep_av<-aggregate(cbind(LooseP,TP_ug)~LakeName,deep_av1,mean,na.rm=TRUE)

# Finally - bind the littoral and averaged deep site dataframe together
# Note: You named another data frame "mix" for the mixed model (Part 4) - so be cognizant of that
mix<-rbind(lit2,deep_av)
tapply(mix$TP_ug,mix$LakeName,mean)
tapply(mix$LooseP,mix$LakeName,mean)


### RAREFACTION ANALYSIS - for TP and loosely-bound P, work by each lake individually.
### SWAN LAKE -------------------------------------------------------------------------------------------------
SW_rare<-subset(mix,LakeName=="Swan")
mean(SW_rare$TP_ug) #1184.913 (whole-lake mean TP - "true" concentration for this analysis)
mean(SW_rare$LooseP) #37.03029 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
sw_sample <- SW_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
sw_mean <- sw_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
sw_mean$RMSE_TP <- sqrt(((sw_mean$mean_TP-1184.913)^2)/1)
mean(sw_mean$RMSE_TP) # add these values to a spreadsheet
sw_mean$RMSE_loose <- sqrt(((sw_mean$mean_loose-37.03029)^2)/1)
mean(sw_mean$RMSE_loose) # add these values to a spreadsheet

### SILVER LAKE -----------------------------------------------------------------------------------------------
SI_rare<-subset(mix,LakeName=="Silver")
mean(SI_rare$TP_ug) #915.5767 (whole-lake mean TP - "true" concentration for this analysis)
mean(SI_rare$LooseP) #167.2061 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
si_sample <- SI_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
si_mean <- si_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
si_mean$RMSE_TP <- sqrt(((si_mean$mean_TP-915.5767)^2)/1)
mean(si_mean$RMSE_TP) # add these values to a spreadsheet
si_mean$RMSE_loose <- sqrt(((si_mean$mean_loose-167.2061)^2)/1)
mean(si_mean$RMSE_loose) # add these values to a spreadsheet

### FIVE ISLAND LAKE ------------------------------------------------------------------------------------------
FI_rare<-subset(mix,LakeName=="Five Island")
mean(FI_rare$TP_ug) #956.6903 (whole-lake mean TP - "true" concentration for this analysis)
mean(FI_rare$LooseP) #108.3087 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
fi_sample <- FI_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
fi_mean <- fi_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
fi_mean$RMSE_TP <- sqrt(((fi_mean$mean_TP-956.6903)^2)/1)
mean(fi_mean$RMSE_TP) # add these values to a spreadsheet
fi_mean$RMSE_loose <- sqrt(((fi_mean$mean_loose-108.3087)^2)/1)
mean(fi_mean$RMSE_loose) # add these values to a spreadsheet

### STORM LAKE ------------------------------------------------------------------------------------------------
SO_rare<-subset(mix,LakeName=="Storm")
mean(SO_rare$TP_ug) #663.1183 (whole-lake mean TP - "true" concentration for this analysis)
mean(SO_rare$LooseP) #47.74041 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
so_sample <- SO_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
so_mean <- so_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
so_mean$RMSE_TP <- sqrt(((so_mean$mean_TP-663.1183)^2)/1)
mean(so_mean$RMSE_TP) # add these values to a spreadsheet
so_mean$RMSE_loose <- sqrt(((so_mean$mean_loose-47.74041)^2)/1)
mean(so_mean$RMSE_loose) # add these values to a spreadsheet

### CENTER LAKE -----------------------------------------------------------------------------------------------
CE_rare<-subset(mix,LakeName=="Center")
mean(CE_rare$TP_ug) #896.9995 (whole-lake mean TP - "true" concentration for this analysis)
mean(CE_rare$LooseP) #71.57667 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
ce_sample <- CE_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
ce_mean <- ce_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
ce_mean$RMSE_TP <- sqrt(((ce_mean$mean_TP-896.9995)^2)/1)
mean(ce_mean$RMSE_TP) # add these values to a spreadsheet
ce_mean$RMSE_loose <- sqrt(((ce_mean$mean_loose-71.57667)^2)/1)
mean(ce_mean$RMSE_loose) # add these values to a spreadsheet

### NORTH TWIN LAKE -------------------------------------------------------------------------------------------
NT_rare<-subset(mix,LakeName=="North Twin")
mean(NT_rare$TP_ug) #1070.393 (whole-lake mean TP - "true" concentration for this analysis)
mean(NT_rare$LooseP) #89.39498 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
nt_sample <- NT_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
nt_mean <- nt_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
nt_mean$RMSE_TP <- sqrt(((nt_mean$mean_TP-1070.393)^2)/1)
mean(nt_mean$RMSE_TP) # add these values to a spreadsheet
nt_mean$RMSE_loose <- sqrt(((nt_mean$mean_loose-89.39498)^2)/1)
mean(nt_mean$RMSE_loose) # add these values to a spreadsheet

### SOUTH TWIN LAKE ---------------------------------------------------------------------------------------
ST_rare<-subset(mix,LakeName=="South Twin")
mean(ST_rare$TP_ug) #935.9872 (whole-lake mean TP - "true" concentration for this analysis)
mean(ST_rare$LooseP) #110.2209 (whole-lake mean loose P - "true" concentration for this analysis)

# Random re-sampling without replacement on a subset of 2 to 9 sampling sites 
# NOTE: keep re-running lines below, but change size of rep_sample_n each time (from 2 to 9)

# random re-sampling, on a subset of n sampling sites, repeated 1000 times
st_sample <- ST_rare %>% 
  rep_sample_n(size = 2, reps = 1000, replace=FALSE)

# calculate the mean concentration (TP and loose P) from each subset (group by replicate, each replicate is a re-sampling)
st_mean <- st_sample %>% 
  ungroup(replicate) %>% 
  mutate(replicate=factor(replicate)) %>% 
  group_by(replicate) %>% 
  summarise(mean_TP = mean(TP_ug),mean_loose = mean(LooseP))

# calculate RMSE for each re-sampling estimate over the 1000 runs, take average
st_mean$RMSE_TP <- sqrt(((st_mean$mean_TP-935.9872)^2)/1)
mean(st_mean$RMSE_TP) # add these values to a spreadsheet
st_mean$RMSE_loose <- sqrt(((st_mean$mean_loose-110.2209)^2)/1)
mean(st_mean$RMSE_loose) # add these values to a spreadsheet

#------------------------------------------------------------------------------------------------------------------
# Read in spreadsheet with mean RMSE values recorded for each subsampling event
rare<-read.csv("Rarefaction_RMSE.csv")

# normalize the RMSE values as a percent of the whole-lake mean P concentration (so that values can be visualized/compared)
rare_norm <- rare %>% 
  mutate(RMSE_norm_TP = (mean_RMSE_TP/mean_TP_whole)*100, RMSE_norm_loose = (mean_RMSE_loose/mean_loose_whole)*100) %>% 
  select(lake,type,n,RMSE_norm_TP,RMSE_norm_loose,mean_RMSE_TP,mean_RMSE_loose)

# compare density plot of raw and normalized RMSE values - make sure distributions look similar, normalizing isn't hiding an unusal pattern
ggplot(data=rare)+
  geom_density(aes(x=mean_RMSE_TP))
ggplot(data=rare_norm)+
  geom_density(aes(x=RMSE_norm_TP))
ggplot(data=rare)+
  geom_density(aes(x=mean_RMSE_loose))
ggplot(data=rare_norm)+
  geom_density(aes(x=RMSE_norm_loose))

# VISUALIZE RESULTS OF RAREFACTION ANALYSIS (Figure 5A-B)
# exclude rows that include all 10 sampling sites (so "true" values, RMSE=0)
rare_sample <- subset(rare_norm, type=="Subsampe")

with(core,levels(LakeName))
brewer.pal(7,"Dark2")
CEcol<-"#A6761D" #brown
FIcol<-"#D95F02" #orange
NTcol<-"#7570B3" #purpleblue
SIcol<-"#E7298A" #pink
SOcol<-"#E6AB02" #mustard
STcol<- "#1B9E77" #teal
SWcol<-"#66A61E" #green

tp_rare<-ggplot(data=rare_sample, aes(x=n, y=RMSE_norm_TP,group=lake,color=lake))+
  geom_point(size=0)+scale_color_manual(values=c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol))+ 
  geom_line(size=1.2)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlab("Number of Sampling Sites") + ylab("Total P Normalized RMSE (%)")+theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=10))+
  annotate("text", x=2.2,y=24,label="A",size=6)+ xlim(2,9) + ylim(0,25) +
  theme(legend.title=element_blank(), legend.position = c(0.7,0.7),legend.text = element_text(size=9),legend.key.size=unit(0.5,"cm"),legend.background=element_blank())

loose_rare<-ggplot(data=rare_sample, aes(x=n, y=RMSE_norm_loose,group=lake,color=lake))+
  geom_point(size=0)+scale_color_manual(values=c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol))+
  geom_line(size=1.2)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlab("Number of Sampling Sites") + ylab("Loosely-Bound P Normalized RMSE (%)")+theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=10))+
  annotate("text", x=2.2,y=24,label="B",size=6)+ xlim(2,9) + ylim(0,25) +
  theme(legend.position = "none")

rare_plots2<-grid.arrange(tp_rare,loose_rare,nrow=1)
ggsave("Figure5.png",rare_plots2,width=5.5,height=3,units="in",dpi=300)


ggplot(data=rare_sample, aes(x=n, y=RMSE_norm_TP,group=lake,color=lake))+
  geom_point(size=0)+scale_color_manual(values=c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol))+ 
  geom_line(size=1.2)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlab("Number of Sampling Sites") + ylab("Total P Normalized RMSE (%)")+theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=10))+
  annotate("text", x=2.2,y=24,label="A",size=6)+ xlim(2,9) + ylim(0,15) +
  theme(legend.title=element_blank(), legend.position = c(0.7,0.7),legend.text = element_text(size=9),legend.key.size=unit(0.5,"cm"),legend.background=element_blank())
ggsave("presentation_rare.png", width=4, height=4, units="in", dpi=300)


ggplot(data=rare_sample, aes(x=n, y=RMSE_norm_TP,group=lake,color=lake))+
  geom_point(size=0)+scale_color_manual(values=c("grey40","#ce1256","grey40","grey40","#ce1256","#ce1256","#ce1256"))+ 
  geom_line(size=1.2)+
  theme_linedraw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  xlab("Number of Sampling Sites") + ylab("Total P Normalized RMSE (%)")+theme(axis.text=element_text(color="black",size=9),axis.title=element_text(size=10))+
  annotate("text", x=2.2,y=24,label="A",size=6)+ xlim(2,9) + ylim(0,15) +
  theme(legend.title=element_blank(), legend.position = "none",legend.text = element_text(size=9),legend.key.size=unit(0.5,"cm"),legend.background=element_blank())
ggsave("presentation_rare2.png", width=4, height=4, units="in", dpi=300)



### PART 7 - VOLUME DEVELOPMENT (FIGURE 4)----------------------------------------------------------------------------------------------------------------------------------------
### PART 7a. Visualize the relationship between the coefficient of variation in loosely-bound P and lake basin volume development (Dv)

# Read in data to calculate Dv
Dv<-read.csv("Dv_Analysis.csv")
View(Dv)

# Unit conversions first. Get lake volumes in cubic meters and surface areas in square meters
Volume_m3<-Dv$Volume_km3*1000000000
SurfaceArea_m2<-Dv$SurfaceArea_ha*10000
# Calculate the volume of a perfect cone (m3). Basal area = lake surface area (m2) and height = lake max depth (m)
ConeVolume<-(0.3333*Dv$Z_max*SurfaceArea_m2)
# Calcualte volume development as the ratio of the lake volume to the volume of the perfect cone with the lake dimensions
Dv_calc<-Volume_m3/ConeVolume
# Bind together columns you'd like for analysis
DV<-cbind(Dv_calc,Dv)

# FIGURE S1
png("FigureS1.png", width=6, height=6, units="in",res=300)
par(mfrow=c(1,1),mar=c(5.1,5.1,4.1,2.1))
plot(DV$CV_Loose~DV$Dv_calc,xlab="Volume Development (Dv)",ylab="Loosely-Bound P Coefficient of Variation (%)",pch=21,cex=3,
     cex.axis=1.2,cex.lab=1.2,xlim=c(0.5,2.5),ylim=c(5,40),las=1,col="black",
     bg=c("goldenrod2","darkslategray4","goldenrod2","goldenrod2","goldenrod2","darkslategray4","darkslategray4"))
text(1.42,40.1,labels="Swan",cex=1.2,adj=0)
text(2,27,labels="North Twin",cex=1.2,adj=0)
text(1.97, 11.5,labels="South Twin",cex=1.2,adj=1)
text(1.95,19,labels="Silver",cex=1.2,adj=1)
text(1.05,37.6,labels="Storm",cex=1.2,adj=0.5)
text(2.09,17.8,labels="Center",cex=1.2,adj=0)
text(0.71,30,labels="Five Island",cex=1.2,adj=0)
dev.off()
