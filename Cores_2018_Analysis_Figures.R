# Project: Variation in Sediment Phosphorus (P) Manuscript - JGR Biogeosciences
# Last Modified 10 December 2020
# Contributers: Ellen Albright, Dr. Grace Wilkinson
# Description: The following script provides all analyses and figures for our manuscript "High inter- and intra-lake variation in sediment phosphorus pools in shallow lakes"

# Data citation: Albright, Ellen; Wilkinson, Grace; Fleck, Rachel; Shingai, Quin (2020): Analysis of Sediment Phosphorus Pools in Shallow Lakes. figshare. Dataset. https://doi.org/10.6084/m9.figshare.13362971.v1 

# The code can be run with the following datasets (NOTE: for each dataset, there is a README file with detailed meta-data)
  # "Cores_ALLDATA_2018.csv" - sediment P chemistry and physical characteristics for 7 shallow lakes in NW Iowa, USA. 
  # "Cores_DepthAgeRelationship.csv" - contains mass accumulation rates and sediment bulk density for study lakes to determine depth-age relationship
  # "Chla_MobileP_Regression.csv" - summarizes mobile (redox-sensitive, loosely-bound, organic) sediment P fractions and long term chlorophyll a concentrations for each lake
  # "Swan_macrophytes_2018.csv" - macrophyte survey across Swan Lake (need "studylakes2018.shp" if you want to map the lake outline)
  # "Swan_SedP_2018.csv" - summary of sediment P results and sampling site coordinates for mapping
  # "Dv_Analysis.csv" - includes morphometric data to calculate lake volume development. Used for Supporting Information Figure.

# Clear environment, set working directory - UPDATE AS NEEDED ----------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
getwd()
setwd('C:/Users/Ellen/Box Sync/DISSERTATION/Albright_Chapter_1/Data/Data_Scripts_SUBMIT')

# REQUIRED PACKAGES FOR CODE - install as needed
# install.packages("tidyverse")
# install.packages("RColorBrewer")
# install.packages("robCompositions")
# install.packages("vegan")
# install.packages("sf")
# install.packages('lmerTest')
# install.packages("ggnewscale")
library(tidyverse) #general
library(RColorBrewer) #general
library(Cairo) #general
library(robCompositions) #PCA
library(vegan) #PCA
library(lmerTest) #mixed model
library(sf) #spatial
library(ggnewscale) #spatial

# COLORBREWER PALLETES, etc.
par(mfrow=c(1,1))
display.brewer.all()

### OVERVIEW OF SCRIPT ORGANIZATION ### -----------------------------------------------------------------------------------------------------------------------------------------------
# PART 1: Descriptions of inter-lake variation in profundal sediment P
#          1a.calculations of various summary statistics (used to make Table 2 in manuscript) 
#          1b.determination of sediment depth-age relationship (used to date core profiles in part 1c.)
#          1c.visualization of total P and P speciation with depth in the sediment profile (Figure 1)
# PART 2: Analysis of inter-lake variation in profundal sediment P
#          2a. compositional data analysis of profundal sediment P speciation among lakes (Figure 2)
# PART 3: Further analysis of water quality implications of inter-lake variation in profundal sediment P
#          3a. simple linear regression of long-term chlorophyll a concentrations as a function of mobile sediment P species (Figure 3)
# PART 4: Intra-lake variation in sediment P - horizontal variation across the lakebed 
#          4a. mixed model regression effects of water depth at the coring location on total sediment P and loosely-bound sediment P by lake (Figure 4)
#          4b. calculations of coefficient of variation for total and loosely-bound P (Table 2)
# PART 5: Macrophyte beds and sediment P Pools
#          5a. Within-lake variation in sediment loosely-bound and total P in Swan Lake (Figure 5)
# PART 6: Supporting information (Figure S1)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### PART 1 - DESCRIPTIONS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (TABLE 2, FIGURE 1) ------------------------------------------------------------------------------------------
# Read in full 2018 sediment dataset and subset out samples from the deep site of the study lakes (SiteType=="Deep")
down<-read.csv("Cores_ALLDATA_2018.csv")
deep<-subset(down,SiteType=="Deep")

### PART 1a. Calculations of summary statistics for Table 2 ----------------------------------------------------------------------------------------------------------------------------
# Mean and standard deviation of total P over all intervals in the deep site core for each lake
tapply(deep$TP_ug,deep$LakeName,mean,na.rm=TRUE)
tapply(deep$TP_ug,deep$LakeName,sd,na.rm=TRUE)

# Average percent contribution of each sediment P fraction
# Need to tidy the dataset a bit before making these calculations because we do not have labile organic P data for all of our lakes
result<-deep %>% 
  filter(SampleID != "CE1B" & SampleID != "CE1F" & SampleID != "NT1A") %>%  #filter out rows that have missing TP data
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,LabileOP,CalciumP,TP_ug) #select relevant columns to work with
result1<-result %>%
  filter(LakeID %in% c("FI","CE","SO","SI")) %>% #first working with 4 lakes that have labile organic P data
  mutate(RefractoryP = TP_ug - (LooseP+IronP+AlP+CalciumP+LabileOP), #calculate refractory organic P concentration and then relative proportions of each fraction
         Per_Loose=LooseP/TP_ug,Per_Fe=IronP/TP_ug,Per_Al=AlP/TP_ug,Per_Ca=CalciumP/TP_ug,Per_Labile=LabileOP/TP_ug,Per_Refrac=RefractoryP/TP_ug,Per_Organic=(RefractoryP+LabileOP)/TP_ug) %>%
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,TP_ug,
         Per_Loose,Per_Fe,Per_Al,Per_Ca,Per_Labile,Per_Refrac,Per_Organic) #add refractory P if you want to calculate percent of total!

tapply(result1$Per_Loose,result1$LakeName,mean) # percent loosely-bound P, averaged across all intervals in the deep site core for each lake
tapply(result1$Per_Fe,result1$LakeName,mean) # percent redox-sensitive P
tapply(result1$Per_Al,result1$LakeName,mean) # percent aluminum-bound P
tapply(result1$Per_Ca,result1$LakeName,mean) # percent calcium-bound P
tapply(result1$Per_Labile,result1$LakeName,mean) # percent labile organic P
tapply(result1$Per_Refrac,result1$LakeName,mean) # percent refractory organic P
tapply(result1$Per_Organic,result1$LakeName,mean) # percent total organic P (as sum of labile and refractory)

result2<-result %>%
  filter(LakeID %in% c("NT","ST","SW")) %>% #next working with 3 lakes that do NOT have labile organic P data
  mutate(OrganicP = TP_ug - (LooseP+IronP+AlP+CalciumP), # total organic P concentration and then relative proportions of each fraction
         Per_Loose=LooseP/TP_ug,Per_Fe=IronP/TP_ug,Per_Al=AlP/TP_ug,Per_Ca=CalciumP/TP_ug,Per_Organic=OrganicP/TP_ug) %>%
  select(LakeName,LakeID,SampleID,TopDepth,BottomDepth,LooseP,IronP,AlP,CalciumP,TP_ug,
         Per_Loose,Per_Fe,Per_Al,Per_Ca,Per_Organic) #add refractory P if you want to calculate percent of total!

tapply(result2$Per_Loose,result2$LakeName,mean) # percent loosely-bound P, averaged across all intervals in the deep site core for each lake
tapply(result2$Per_Fe,result2$LakeName,mean) # percent redox-sensitive P
tapply(result2$Per_Al,result2$LakeName,mean) # percent aluminum-bound P
tapply(result2$Per_Ca,result2$LakeName,mean) # percent calcium-bound P
tapply(result2$Per_Organic,result2$LakeName,mean) # percent total organic P (assumed to be difference between total P and 4 measured fractions)

# Calculate percent change in sediment total P over time (i.e. over the sediment profile)
# % Increase=new number (surface, 1A) - old (bottom, 1G). Divide increase by original number.
# SWAN: 91.5% increase over time (nearly doubled)
((deep[deep$LakeName=="Swan"&deep$SampleID=="SW1A","TP_ug"]-deep[deep$LakeName=="Swan"&deep$SampleID=="SW1G","TP_ug"]))/deep[deep$LakeName=="Swan"&deep$SampleID=="SW1G","TP_ug"]
# CENTER 66.32%
(deep[deep$LakeName=="Center"&deep$SampleID=="CE1A","TP_ug"]-deep[deep$LakeName=="Center"&deep$SampleID=="CE1G","TP_ug"])/deep[deep$LakeName=="Center"&deep$SampleID=="CE1G","TP_ug"]
# NORTH TWIN 137.1% (more than doubled!)
(deep[deep$LakeName=="North Twin"&deep$SampleID=="NT1B","TP_ug"]-deep[deep$LakeName=="North Twin"&deep$SampleID=="NT1G","TP_ug"])/deep[deep$LakeName=="North Twin"&deep$SampleID=="NT1G","TP_ug"]
# FIVE ISLAND 38.4%
(deep[deep$LakeName=="Five Island"&deep$SampleID=="FI1A","TP_ug"]-deep[deep$LakeName=="Five Island"&deep$SampleID=="FI1G","TP_ug"])/deep[deep$LakeName=="Five Island"&deep$SampleID=="FI1G","TP_ug"]
# SILVER 21.7%
(deep[deep$LakeName=="Silver"&deep$SampleID=="SI1A","TP_ug"]-deep[deep$LakeName=="Silver"&deep$SampleID=="SI1G","TP_ug"])/deep[deep$LakeName=="Silver"&deep$SampleID=="SI1G","TP_ug"]
# STORM 0.8%
(deep[deep$LakeName=="Storm"&deep$SampleID=="SO1A","TP_ug"]-deep[deep$LakeName=="Storm"&deep$SampleID=="SO1G","TP_ug"])/deep[deep$LakeName=="Storm"&deep$SampleID=="SO1G","TP_ug"]
# SOUTH TWIN 17.4%
(deep[deep$LakeName=="South Twin"&deep$SampleID=="ST1A","TP_ug"]-deep[deep$LakeName=="South Twin"&deep$SampleID=="ST1G","TP_ug"])/deep[deep$LakeName=="South Twin"&deep$SampleID=="ST1G","TP_ug"]

### PART 1b. Determination of sediment depth-age relationship (used to date core profiles in Figure 1) ---------------------------------------------------------------------------------
age<-read.csv("Cores_DepthAgeRelationship.csv") 
   # Contains modern (post 1850) sediment mass accumulation rates (MAR) in g/cm2/year from Heathcote et al. (2013) 
   # Also added information on depth intervals of interest for each study lake included in Heathcote's work, plus dry bulk density (BD, g/cm3) for each slice (from phys dataset, 2018)
# STEP 1: Calculate linear sedimentation rate (cm/year): Divide MAR/BD. g/cm2/year * cm3/g = cm/year
LSR<-age$Modern_MAR/age$BD # tells you how many cm of sediment were laid down each year for a given depth interval
# STEP 2: Use LSR to determine the age-range of each core slice (if x cm of sediment were laid down a year, how many years do these 2cm represent?)
Age_Range<-age$Interval/LSR
View(age)
# STEP 3: Bind it all together and calculate the end date of the first slice
age2<-cbind(age,LSR,Age_Range)
age2$End<-age2$Start-age2$Age_Range #keep re-running this line after you manually update start cell values (below)
age2[2,7]<-2001.862 
age2[3,7]<-1986.000
age2[4,7]<-1969.600
age2[5,7]<-1952.877	
age2[6,7]<-1936.062
age2[7,7]<-1893.370
age2[9,7]<-1997.283
age2[10,7]<-1976.094
age2[11,7]<-1954.471
age2[12,7]<-1932.867
age2[13,7]<-1911.093
age2[14,7]<-1855.999
age2[16,7]<-2008.783
age2[17,7]<-1999.392
age2[18,7]<-1989.827
age2[19,7]<-1980.262
age2[20,7]<-1970.523
age2[21,7]<-1945.958
age2[23,7]<-2003.093
age2[24,7]<-1987.553
age2[25,7]<-1972.323
age2[26,7]<-1956.733
age2[27,7]<-1941.143
age2[28,7]<-1902.603
age2[29,7]<-2018
age2[30,7]<-2000.293
age2[31,7]<-1982.586
age2[32,7]<-1964.569
age2[33,7]<-1946.672
age2[34,7]<-1928.551
age2[35,7]<-1882.344
age2[37,7]<-2003.467
age2[38,7]<-1988.505
age2[39,7]<-1973.491
age2[40,7]<-1958.498
age2[41,7]<-1943.427
age2[42,7]<-1904.649
age2[44,7]<-2003.623
age2[45,7]<-1988.784
age2[46,7]<-1973.730
age2[47,7]<-1958.363
age2[48,7]<-1942.755
age2[49,7]<-1903.932
age2$End<-age2$Start-age2$Age_Range

### PART 1c. Visualization of total P and P speciation with depth in the sediment profile (Figure 1)-----------------------------------------------------------------------------------
# This figure illustrates the concentrations of loosely-bound, redox-sensitive, labile organic, Al-bound, Ca-bound, and total P with depth in the sediments
# Colors used for each fraction
cal = "#51AAAE"
al = "#3288bd"
loose = "#d53e4f"
labile = "#fdae61"
iron = "#f46d43"

Cairo(file="Figure1.png", 
      type="png",
      units="in", 
      width=7.5, 
      height=10, 
      pointsize=12, 
      dpi=300)
dev.off()

#Silver Lake
#windows(height = 9, width = 6.5)
par(mfrow=c(4,2), omi=c(0.3,0.3,0.3,0.3), mai = c(0.2,0.3,0.1,0.3))
plot((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","CalciumP"]),
     type="o", pch=21, bg="white",col=cal, cex=1.5, lwd=2, 
     ylim=c(20,0), xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="", ylab = "", las=2)
axis(side=3, at=c(0,400,800,1200,1600), labels=c("0","400","800","1200","1600"), cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),labels=c("2018","1997","1976","1954","1933","1911","1856","1800"), cex.axis=1.2, las=2)
points((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","LabileOP"]),
       type="o",pch=23,bg="white",col=labile,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Silver","BottomDepth"])~(deep[deep$LakeName=="Silver","TP_ug"]),
       type="o",pch=20,col="black",cex=2,lwd=2)
text(1400,0.5,labels="Silver",cex=1.5,adj=c(0,NULL))

#Center Lake
plot((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2,xaxt='n',xlab="",yaxt='n', las=2)
axis(side=3, at=c(0,400,800,1200,1600),labels=c("0","400","800","1200","1600"),cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2002","1986","1970","1953","1936","1893","1850"),
     cex.axis=1.2, las=2)
points((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","LabileOP"]),
       type="o",pch=23,bg="white",col=labile,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Center","BottomDepth"])~(deep[deep$LakeName=="Center","TP_ug"]),
       type="o",pch=20,col="black",cex=2,lwd=2)
text(1370,0.5,labels="Center",cex=1.5,adj=c(0,NULL))

#Five Island
plot((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="", las=2)
#axis(side=3, at=c(0,400,800,1200,1600), labels=c("0","400","800","1200","1600"),cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2009","1999","1990","1980","1971","1946","1921"),
     cex.axis=1.2, las = 2)
points((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","LabileOP"]),
       type="o",pch=23,bg="white",col=labile,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Five Island","BottomDepth"])~(deep[deep$LakeName=="Five Island","TP_ug"]),
       type="o",pch=20,col="black",cex=2,lwd=2)
text(1170,0.5,labels="Five Island",cex=1.5,adj=c(0,NULL))

#Storm
plot((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="", yaxt="n")
#axis(side=3, at=c(0,400,800,1200,1600),labels=c("0","400","800","1200","1600"),cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2003","1988","1972","1957","1941","1903","1864"),
     cex.axis=1.2, las=2)
points((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","LabileOP"]),
       type="o",pch=23,bg="white",col=labile,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Storm","BottomDepth"])~(deep[deep$LakeName=="Storm","TP_ug"]),
       type="o",pch=20,col="black",cex=1.8,lwd=2)
text(1400,0.5,labels="Storm",cex=1.5,adj=c(0,NULL))

#North Twin
plot((deep[deep$LakeName=="North Twin","BottomDepth"])~(deep[deep$LakeName=="North Twin","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="", las=2)
#axis(side=3, at=c(0,400,800,1200,1600), 
     #labels=c("0","400","800","1200","1600"),
     #cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2000","1982","1964","1946","1928","1882","1835"),
     cex.axis=1.2, las=2)
points((deep[deep$LakeName=="North Twin","BottomDepth"])~(deep[deep$LakeName=="North Twin","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="North Twin","BottomDepth"])~(deep[deep$LakeName=="North Twin","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="North Twin","BottomDepth"])~(deep[deep$LakeName=="North Twin","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="North Twin","BottomDepth"])~(deep[deep$LakeName=="North Twin","TP_ug"]),
       type="o",pch=20,col="black",cex=1.8,lwd=2)
text(1180,0.5,labels="North Twin",cex=1.5,adj=c(0,NULL))

#South TWin
plot((deep[deep$LakeName=="South Twin","BottomDepth"])~(deep[deep$LakeName=="South Twin","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="",yaxt='n')
#axis(side=3, at=c(0,400,800,1200,1600), 
     #labels=c("0","400","800","1200","1600"),cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2003","1988","1973","1958","1943","1904","1865"), col.axis="gray50",cex.axis=1.2, las=2)
points((deep[deep$LakeName=="South Twin","BottomDepth"])~(deep[deep$LakeName=="South Twin","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="South Twin","BottomDepth"])~(deep[deep$LakeName=="South Twin","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="South Twin","BottomDepth"])~(deep[deep$LakeName=="South Twin","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="South Twin","BottomDepth"])~(deep[deep$LakeName=="South Twin","TP_ug"]),
       type="o",pch=20,col="black",cex=1.8,lwd=2)
text(1170,0.5,labels="South Twin",cex=1.5,adj=c(0,NULL))

#Swan 
plot((deep[deep$LakeName=="Swan","BottomDepth"])~(deep[deep$LakeName=="Swan","CalciumP"]),
     type="o",pch=21,bg="white",col=cal,cex=1.5,lwd=2, 
     ylim=c(20,0),xlim=c(1,1650), cex.axis=1.2, xaxt='n', xlab="",las=2)
#axis(side=3, at=c(0,400,800,1200,1600), 
     #labels=c("0","400","800","1200","1600"),cex.axis=1.2)
axis(side=4,at=c(0,2,4,6,8,10,15,20),
     labels=c("2018","2003","1989","1974","1958","1943","1904","1865"), col.axis="gray50",cex.axis=1.2, las=2)
points((deep[deep$LakeName=="Swan","BottomDepth"])~(deep[deep$LakeName=="Swan","AlP"]), 
       type="o", pch=25,bg="white",col=al, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Swan","BottomDepth"])~(deep[deep$LakeName=="Swan","LooseP"]),
       type="o",pch=24,bg="white",col=loose, cex=1.5, lwd=2)
points((deep[deep$LakeName=="Swan","BottomDepth"])~(deep[deep$LakeName=="Swan","IronP"]),
       type="o",pch=22,bg="white",col=iron,cex=1.5,lwd=2)
points((deep[deep$LakeName=="Swan","BottomDepth"])~(deep[deep$LakeName=="Swan","TP_ug"]),
       type="o",pch=20,col="black",cex=1.8,lwd=2)
text(1400,0.5,labels="Swan",cex=1.5,adj=c(0,NULL))

plot(1,1,col='white',xaxt='n',yaxt='n',bty='n')
legend("center",legend=c('Aluminum-Bound','Calcium-Bound','Loosely-bound', 'Redox-senitive','Labile Organic','Total'), pch=c(25,21,24,22,23,16),col=c(al, cal, loose, iron, labile, "black"),pt.bg=c("white","white","white","white","white","black"),lwd=3, pt.cex=2, cex=1.4, y.intersp=1,bty='n')


### PART 2 - ANALYSIS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (TABLE 2, FIGURE 1) ---------------------------------------------------------------------------------------------
### PART 2a. compositional data analysis (CoDA) of profudnal sediment P speciation among lakes (Figure 2) -----------------------------------------------------------------------------
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
# Take the centered log raio (CLR) of the composiitonal data
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

# PCA Bioplot (Figure 2)
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

par(mfrow=c(1,1),mar=c(4,4,0.8,0.8)) #bottom, left, top, right
biplot(core.rda,choices=c(1,2),display="sites",type="points",scaling=3,xlab="",ylab="",ylim=c(-0.8,0.8),xlim=c(-0.8,0.8)) #plot PCA, symmetrical scaling, "sites" only
title(xlab="PC1 (45.27)%",ylab="PC2 (36.42%)",cex.lab=1.2,) #show % variation explained by each principal component
with(core,points(core.rda,display="sites",col="black",pch=21,scaling=3,bg=colvec[LakeName],cex=1.5)) #color points by lake
ordihull(core.rda,group=core$LakeName,scaling=3,label=FALSE,draw="polygon",col=c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol),alpha=0.2)
text(-0.7,-0.31,labels="Swan",col="#3e6512",cex=1.2)
text(0.2,-0.7,labels="North Twin",col="#4f4a8c",cex=1.2)
text(0.45, -0.5,labels="South Twin",col="#105d46",cex=1.2)
text(0.5,0.1,labels="Silver",col="#b01463",cex=1.2)
text(0.13,0.6,labels="Storm",col="#9a7301",cex=1.2)
text(-0.2,0.3,labels="Center",col="#654812",cex=1.2)
text(-0.2,-0.08,labels="Five Island",col="#8d3e01",cex=1.2)

# envfit() function to add environmental data 
x_clr<-envfit(core.rda,core_CLR[,12:16],choices=c(1,2))
plot(x_clr,choices=c(1,2),at=c(0,0),axis=FALSE,p.max=NULL,add=TRUE,col="black",labels=c(" "," "," "," "," "),cex=1.5)
text(-0.5,0.32,labels="Redox-
     Sensitive P")
text(-0.38,0.07,labels="Al-Bound P")
text(-0.21,-0.5,labels="Organic P")
text(0.65,-0.41,labels="Loosely-Bound P")
text(0.67,0.53,labels="Ca-Bound P")

Cairo(file="Figure2.png", 
      type="png",
      units="in", 
      width=7, 
      height=7, 
      pointsize=12, 
      dpi=300)
dev.off()

### PART 3 - FURTHER ANALYSIS OF WATER QUALITY IMPLICATIONS OF INTER-LAKE VARIATION IN PROFUNDAL SEDIMENT P (FIGURE 3) ----------------------------------------------------------------
### PART 3a. simple linear regression of long-term chlorophyll a concentrations as a function of mobile sediment P species (Figure 3) -------------------------------------------------
# Read in auxiliary dataset that summarizes mobile sediment P species and long-term chla concentrations
chla<-read.csv("Chla_MobileP_Regression.csv")
View(chla)

x<-chla$mobilePercent*100
SDx<-(chla$mobileSD/sqrt(chla$mobileN))*100
xlabel<-"% Mobile P of Total Sediment P"

y<-chla$chlorophyll
SDy<-chla$chlorophyllSD/sqrt(chla$chlorophyllN)

windows(height = 4, width = 4)
par(mai=c(1,1,0.1,0.1))
plot(x,y, pch = 19, cex=1.5, xlim = c(80,100), ylim = c(0,100),
     xlab = xlabel, ylab = expression(Average~Chlorophyll~"("*mu*g~L^-1*")"), cex.lab=1)
arrows(x, y-SDy, x1=x, y1=y+SDy, angle=0, length=0)
arrows(x-SDx, y, x1=x+SDx, y1=y, angle=0, length=0)
mod<-lm(y~x)
summary(mod) #Adj R2=0.5232, F(1,5)=7.584, p=0.04013
lines(c(82,96), c(82*3.691-275.227, 96*3.691-275.227), col="dodgerblue3", lwd=3)

Cairo(file="Figure3.png", 
      type="png",
      units="in", 
      width=7, 
      height=7, 
      pointsize=12, 
      dpi=300)
dev.off()

### PART 4 - INTRA-LAKE VARIATION IN SEDIMENT P - HORIZONTAL VARIATION ACROSS THE LAKEBED (FIGURE 4, TABLE 2) -------------------------------------------------------------------------
### PART 4a. Mixed model regression effects of water depth at the coring location on total sediment P and loosely-bound sediment P by lake (Figure 4) ---------------------------------

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
#Calculate the variance explained by the model (R-squared)
rsq(mix$TP_ug, predict(totalP_mixed)) #R2=0.631
# Null model - don't allow lake to be a random effect on the intercept
totalP_null<-lm(TP_ug ~ WaterDepth,data=mix)
summary(totalP_null)
#Compare the mixed model to the null model
anova(totalP_mixed, totalP_null) #including random effects on model intercept improved fit (p<0.0001)

# PLOT THE DATA ============================================================
windows(height=6.5, width=4)
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

Cairo(file="Figure4.png", 
      type="png",
      units="in", 
      width=4, 
      height=7, 
      pointsize=12, 
      dpi=300)
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

### PART 5 - INTRA-LAKE VARIATION IN SEDIMENT P - MACROPHYTE BEDS AND SEDIMENT P POOLS Macrophyte beds and sediment P Pools (Figure 5) ------------------------------------------------
### PART 5a. Within-lake variation in sediment loosely-bound and total P in Swan Lake (Figure 5)---------------------------------------------------------------------------------------

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

# Plot data together (Figure 5)
ggplot() + 
  geom_sf(data = swan_shp, size = 1, color = "black", fill = "white") + #f2f8f9 #e0eff1
  geom_sf(data=sites,aes(size=TOTAL),col="#b8e186")+scale_size(range=c(5,11))+
  new_scale("size")+
  geom_sf(data=cores,aes(size=TP_ug,fill=LoosePER),shape=21,col="black")+scale_size(range=c(5,16))+scale_fill_continuous(low=("#fde0ef"),high=("#8e0152"))+
  theme_classic()+
  coord_sf()

Cairo(file="Figure5.png", 
      type="png",
      units="in", 
      width=8, 
      height=7, 
      pointsize=12, 
      dpi=300)
dev.off()

### PART 6 - SUPPORTING INFORMATION (FIGURE S1)----------------------------------------------------------------------------------------------------------------------------------------
### PART 6a. Visualize the relationship between the coefficient of variation in loosely-bound P and lake basin volume development (Dv)

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
par(mfrow=c(1,1),mar=c(5.1,5.1,4.1,2.1))
plot(DV$CV_Loose~DV$Dv_calc,xlab="Volume Development (Dv)",ylab="Loosely-Bound P Coefficient of Variation (%)",pch=21,cex=3,
     cex.axis=1.2,cex.lab=1.2,xlim=c(0.5,2.5),ylim=c(5,40),las=1,col="black",
     bg=c("goldenrod2","darkslategray4","goldenrod2","goldenrod2","goldenrod2","darkslategray4","darkslategray4"))
text(1.42,40.1,labels="Swan",cex=1.2,adj=0)
text(2.19,24.3,labels="North Twin",cex=1.2,adj=0)
text(1.97, 11.5,labels="South Twin",cex=1.2,adj=1)
text(1.99,20.76,labels="Silver",cex=1.2,adj=1)
text(1.05,37.6,labels="Storm",cex=1.2,adj=0.5)
text(2.09,17.8,labels="Center",cex=1.2,adj=0)
text(0.71,29.3,labels="Five Island",cex=1.2,adj=0)

Cairo(file="FigureS1.png", 
      type="png",
      units="in", 
      width=7, 
      height=7, 
      pointsize=12, 
      dpi=300)
dev.off()