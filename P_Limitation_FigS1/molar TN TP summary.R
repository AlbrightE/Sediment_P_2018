# Analysis of the molar ratio of TN:TP
# 7 Iowa Lakes for Albright et al. 

### NOTE: csv is named "Study Lakes 2000-2020.csv" instead of "Iowa Lakes Dataset 2000-2020.csv" Sorry.

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)
if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

# CENTER ================================================
center = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='Center Lake') %>%
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  # Let's do some combining of data in columns to fix this unholy mess
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ 
                           `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                        TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(center$molarNP, na.rm = TRUE)

# FIVE ISLAND ================================================
five = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='Five Island Lake') %>%
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  # Let's do some combining of data in columns to fix this unholy mess
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                        TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(five$molarNP, na.rm = TRUE)

# NORTH TWIN ================================================
ntwin = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='North Twin Lake (maximum water depth)') %>%
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                        TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(ntwin$molarNP, na.rm = TRUE)

# SILVER LAKE ================================================
silver = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='Silver Lake Max Depth') %>%
  filter(county=="Dickinson") %>%
  select(-facilityID:-name, -huc8:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                        TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(silver$molarNP, na.rm = TRUE)


# STORM LAKE ================================================
storm = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='Storm Lake Max Depth') %>%
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  # Let's do some combining of data in columns to fix this unholy mess
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                        TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(storm$molarNP, na.rm = TRUE)

# SWAN LAKE ================================================
swan = read.csv("Iowa Lakes Dataset 2000-2020.csv") %>%
  filter(name=='Swan Lake Max Depth') %>%
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -unit:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  mutate(TKN = case_when(is.na(.$`Nutrient-nitrogen`) ~ `Kjeldahl nitrogen (AS N)`, 
                         TRUE ~ `Nutrient-nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  mutate(TP = case_when(is.na(.$`Phosphate-phosphorus (as P)`) ~ 
                          `Phosphate-phosphorus (AS P)`, 
                         TRUE ~ `Phosphate-phosphorus (as P)`)) %>%
  filter(TKN>0) %>%
  select(-'Depth, bottom':-Phycocyanin) %>%
  mutate(TN = TKN + NOx,
         molarNP = ((TN/1000)/14.01) / ((TP/1000)/30.97))
mean(swan$molarNP, na.rm = TRUE)

# SOUTH TWIN ==============================================
stwin = read.csv("South Twin Data.csv") %>%
  # Let's clean this mess up and pivot
  select(-facilityID:-sampleDate, -projectCode, -cas_rn, 
         -fraction, -relDepth, -quantFlag:-remark) %>%
  pivot_wider(id_cols = c(year, doy), names_from = analyte, values_from = result) %>%
  filter(year<2018) %>%
  # Let's do some combining of data in columns to fix this unholy mess
  mutate(TKN = case_when(is.na(.$`Kjeldahl nitrogen`) ~ 
                                `Kjeldahl nitrogen (AS N)`, 
                              TRUE ~ `Kjeldahl nitrogen`)) %>%
  mutate(NOx = case_when(is.na(.$`Inorganic nitrogen (nitrate and nitrite) (as N)`) ~ 
                           `Inorganic nitrogen (nitrate and nitrite) (AS N)`, 
                         TRUE ~ `Inorganic nitrogen (nitrate and nitrite) (as N)`)) %>%
  rename(TP = "Phosphate-phosphorus (AS P)") %>%
  select(-`Ammonia-nitrogen (as N)`:-`Inorganic nitrogen (nitrate and nitrite) (as N)`,
         -`Inorganic nitrogen (nitrate and nitrite) (AS N)`:-`Kjeldahl nitrogen`) %>%
  #Now we can finally calculate molar ratios - the science. Amen.
  mutate(molarNP = (((TKN + NOx)/1000)/14.01) / ((TP/1000)/30.97))
mean(stwin$molarNP, na.rm = TRUE)

# Colors
brewer.pal(7,"Dark2")
CEcol<-"#A6761D" #brown
FIcol<-"#D95F02" #orange
NTcol<-"#7570B3" #purpleblue
SIcol<-"#E7298A" #pink
SOcol<-"#E6AB02" #mustard
STcol<- "#1B9E77" #teal
SWcol<-"#66A61E" #green
colvec<-c(CEcol,FIcol,NTcol,SIcol,STcol,SOcol,SWcol)

windows(height = 4, width = 6.5)
par(mai = c(1,1,0.2,0.2))
boxplot(center$molarNP, five$molarNP, ntwin$molarNP,
        silver$molarNP, stwin$molarNP, storm$molarNP, swan$molarNP,
        col = colvec, ylim = c(0,400), pch = 20,
        names = c("Center", "Five\nIsland", "N Twin", "Silver",
                  "S Twin", "Storm", "Swan"), las = 2, lty = 1,
        ylab = "Molar TN:TP", cex.lab = 1.5)
abline(16, 0, lty = 2, lwd = 2)

