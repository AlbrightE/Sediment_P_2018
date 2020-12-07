# Project: Variation in Sediment Phosphorus (P) Manuscript - JGR Biogeosciences
# Last Modified: 7 December 2020 Ellen Albright
# Description: The following code outlines our calculations and data claning procedure to determine sediment physical characteristics and P concentrations

# This script can be run with the following example datasets:
  # Sediment physical data: "Cores_Physical_EXAMPLE.csv" (meta-data in "Cores_Physical_EX_README.xlsx")
  # Sediment P speciation: "Cores_SedP_EXAMPLE.csv" (meta-data in "Cores_SedP_EX_README.xlsx")
  # Sediment Total P: "Cores_TP_EXAMPLE.csv" (meta-data in "Cores_TP_EX_README.xlsx")

# Libraries needed:
# install.packages("tidyverse")
library(tidyverse)

# Set working directory - UPDATE AS NEEDED
getwd()
setwd('C:/Users/Ellen/Box Sync/DISSERTATION/Albright_Chapter_1/Data/Data_Scripts_SUBMIT')

##### STEP 1 - Calculations for SEDIMENT PHYSICAL CHARACTERISTICS (moisture content, organic matter, bulk density) -------------------------------------------------------------------------------------------
# Read in example dataset
rawphys<-read.csv('Cores_Physical_EXAMPLE.csv',as.is=T)

# Calculations
# 1a. Dry Mass (g)
dry<-rawphys$DryMassTin-rawphys$TinTare
# 1b. Moisture Content (fraction) ***These values are needed for later to calculate the dry mass equivalent of the fresh sediment used in the sequential P extractions (See Step 2 Section)
MCfrac<-(rawphys$physWetMass-(rawphys$DryMassTin-rawphys$TinTare))/rawphys$physWetMass
# 1c. Combusted Mass (g)
combusted<-(rawphys$CombustedMassTin-rawphys$TinTare)
# 1d. Organic Matter Content (%) as loss on ignition (LOI)
OM<-(((rawphys$DryMassTin-rawphys$TinTare)-(rawphys$CombustedMassTin-rawphys$TinTare))/(rawphys$DryMassTin-rawphys$TinTare))*100
# 1e. Bulk Density (g/cm3) - equation from Hakanson & Jansson 2002
bd<-260/(100+(1.6*((MCfrac*100)+(OM/100)*(100-(MCfrac*100)))))

# Bind all columns to write complete dataset of sediment physical characteristics ('phys')
phys<-cbind(rawphys,dry,MCfrac,combusted,OM,bd)

##### STEP 2 - Calculations for CONCENTRATIONS OF SEDIMENT P FRACTIONS (loosely-bound, redox-sensitive, labile organic, calcium-bound, aluminum-bound) ---------------------------------------------------------
# Read in example dataset
rawsed<-read.csv("Cores_SedP_EXAMPLE.csv",as.is=T)

# 2a. Calculate the DRY MASS EQUIVALENT of the fresh sediment subsample used for the sequential P extractions (units=g)
DryMassEq<-rawsed$WetMass*(1-rawsed$Mcfrac)

# 2b. Calculate the concentration of LOOSELY-BOUND P, NH4Cl Extraction
NH4Cl_SRP_Corrected<-(rawsed$NH4Cl_SRP1+rawsed$NH4Cl_SRP2)/2 #average of SRP concentration of spectrophotometer duplicates (units=ug/L)
LooseP<-(NH4Cl_SRP_Corrected*rawsed$VolNH4Cl)/DryMassEq #concentration of Loosely-Bound P (units=ugP/g dry sediment)

# 2c. Calculate the concentration of REDOX-SENSITIVE P, BD Extraction 
BD_SRP_Corrected<-(rawsed$BD_SRP1+rawsed$BD_SRP2)/2 #average of SRP concentration of spectrophotometer duplicates (units=ug/L)
IronP<-(BD_SRP_Corrected*rawsed$VolBD)/DryMassEq #concentration of Loosely-Bound P (units=ugP/g dry sediment)

# 2d. Calculate the concentration of ALUMINUM-BOUND & LABILE ORGANIC P, NaOH Extraction 
NaOH_SRP_Corrected<-(rawsed$NaOH_SRP1+rawsed$NaOH_SRP2)/2 #average SRP concentration of NON-DIGESTED spectrophotometer duplicates (units=ug/L)
NaOH_pHSRP<-(NaOH_SRP_Corrected*(rawsed$NaOH_postPH-rawsed$NaOH_JarTare)/(rawsed$NaOH_prePH-rawsed$NaOH_JarTare))/1.00152 #correct the average SRP concentration for pH adjustment (units=ug/L)
AlP<-(NaOH_pHSRP*rawsed$VolNaOH)/DryMassEq #Concentration of Al-Bound P (units=ugP/g dry sediment)

NaOH_digSRP_Corrected<-(rawsed$NaOH_digSRP1+rawsed$NaOH_digSRP2)/2 #average SRP concentration of DIGESTED spectrophotometer duplicates (units=ug/L)
NaOH_digpHSRP<-(NaOH_digSRP_Corrected*(rawsed$NaOH_postPH-rawsed$NaOH_JarTare)/(rawsed$NaOH_prePH-rawsed$NaOH_JarTare))/1.00152 #correct the average SRP concentration for pH adjustment(units=ug/L)
NaOH_extractableP<-(NaOH_digpHSRP*rawsed$VolNaOH)/DryMassEq #Concentration of all NaOH-extractable P (units=ugP/g dry sediment)
LabileOP<-NaOH_extractableP-AlP #Concentration of labile organic P (units=ugP/g dry sediment)

# 2e. Calcuate the concentration of CALCIUM-BOUND P, HCl Extraction 
HCl_SRP_Corrected<-(rawsed$HCl_SRP1+rawsed$HCl_SRP2)/2 #average of SRP concentration of spectrophotometer duplicates (units=ug/L)
HCl_pHSRP<-(HCl_SRP_Corrected*((rawsed$HCl_postPH-rawsed$HCl_JarTare)/(rawsed$HCl_prePH-rawsed$HCl_JarTare)))/1.00452 #pH-corrected, average SRP concentration (ug/L)
CalciumP<-(HCl_pHSRP*rawsed$VolHCl)/DryMassEq #Concentration of Calcium-Bound P (ugP/g dry sediment)

# Bind all columns to write complete dataset of sediment P fractions ('sed_frac')
sed_frac<-cbind(rawsed,DryMassEq,LooseP,IronP,AlP,LabileOP,CalciumP)
names(sed_frac)

##### STEP 3 - Caclulations for SEDIMENT TOTAL P concentration ------------------------------------------------------------------------------------------------------------------------------------------------
# Read in example dataset
rawTP<-read.csv("Cores_TP_EXAMPLE.csv",as.is=T)

# Calculations
# 3a. Average TP concentration (average dups), corrected for dilution factor (multiply by) (units=mg/L)
avTPdil<-((rawTP$TP1+rawTP$TP2)/2)*rawTP$TP_DilFac
# 3b. Correct conentration from 3a. for pH adjustment (units=mg/L)
CorTP<-(avTPdil*((rawTP$TP_postPH-rawTP$TP_Beaker)/(rawTP$TP_prePH-rawTP$TP_JarTare)))/1.00452
# 3c. Calculate sediment TP concentration (units=mg P/g dry sediment)
TP_mg<-(CorTP*rawTP$TP_DilVol)/rawTP$DryMass
# 3d. Convert calculation from 3c. to ug P/g dry sediment
TP_ug<-TP_mg*1000

#Bind it all together
tp<-cbind(rawTP,TP_mg,TP_ug)

##### STEP 4 - Write your FINAL DATAFRAME based on the desired variables --------------------------------------------------------------------------------------------------------------------------------------

# Select Columns that you want (you'll need dplyr for this!)
# Make sure you select a common variable to join dataframes later. I used "SampleID"
phys2<-select(.data=phys,LakeName,SiteID,LakeID,SampleID,SiteType,WaterDepth,TopDepth,BottomDepth,MCfrac,OM,bd) 
sed2<-select(.data=sed_frac,SampleID,DryMassEq,LooseP,IronP,AlP,LabileOP,CalciumP)
TP2<-select(.data=tp,SampleID,TP_ug)

data<-full_join(phys2,sed2,by="SampleID")
data2<-full_join(data,TP2,by="SampleID")
View(data2) #make sure everything looks good!

write.csv(data2,file="Cores_ALLDATA_2018.csv")



#References
#1. Hakanson, L., & Jansson, M. (2002). Principles of lake sedimentology. Caldwell, NJ: The Blackburn Press. https://doi.org/10.1002/iroh.19850700318 