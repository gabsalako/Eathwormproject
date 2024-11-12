library(raster) #Work with raster data
library(rgdal) #Export GeoTIFFs and other core GIS functions
library (randomForest)
library(randomForestSRC)
library(rpart)
library(sp)	
library(sf)
library(dplyr)
library(caret)#for confusion matrix
library(corrplot)	
library(cluster)
library(caTools)#for data split
library(biomod2)
library(caret)#for confusion matrix
library(usdm)#for VIF
library(spatialRF)
library(pdp)
library(dismo)#for GBM dismo
library(mapview) # for view data on openstreet and other mapping platforms
library(ggpubr)# for ggpar
library(gbm)#for variable contribution, gbm plot
library(gstat)
library(ecospat) # to plot corr matrix with graph
library(mgcv) #to run GAM
library(moments)
library(dplyr)
library(GGally)# for ggpair, rescale01
library(rms) # for C AUC
library(ggplot2)
library(pROC)#for evaluating model e.g AUC(ROC)'
library(SSDM)
library(rms) # for C AUC
library(car)#for r2 scatter plot alternative to ggplot r2
library(lattice)#also for ggplot2
library(mapview)
library(reshape2)#for melt function
library(plyr)#mean sd. se
library(goeveg)#species response curve
library(gganimate)

library(terra)
library(tidyr)
library(timetk)

setwd("C:/Users/gsalako/Documents/BonaRes_Project/DynamicSDM_R")
list.files()


proj4Str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
EWsdym <- read.csv("SDMdy_EndogDens.csv")
head (EWsdym, 20)

EWsdym[,'Observation.date..Sampling.event.'] = as.factor(EWsdym[,'Observation.date..Sampling.event.'])#categorical to display indiv. date
ggplot(data = EWsdym, aes(x = Observation.date..Sampling.event., y =  Density, group=1)) + # if only one data per date
  geom_line() +
  geom_smooth(se = TRUE) +
  labs(x = "Year", y = "Abundance") 

p <- aggregate(EWsdym$Density, by=list(EWsdym$Observation.date..Sampling.event.), FUN=mean)
head(p, 20)
p[,'Group.1'] = as.factor(p[,'Group.1'])#categorical

ggplot(data = p, aes(x =Group.1, y = x,group=1)) + # if only one data per date
  geom_line() +
  geom_smooth(se = TRUE) +
  labs(x = "Year", y = "Abundance") 


setwd("C:/Users/gsalako/Documents/BonaRes_Project/DynamicSDM_R")
list.files()


Endgsdym <- read.csv("Abg_SDMdyEndog.csv" )# for line trends
head(Endgsdym)
str(Endgsdym)
DFEndgDy <- data.frame(Endgsdym)
head(DFEndgDy)
Endgsdym[,'Observation.date..Sampling.event.'] = as.factor(Endgsdym[,'Observation.date..Sampling.event.'])#categorical to display indiv. date
P1 <-ggplot(data = Endgsdym, aes(x =Observation.date..Sampling.event., y = Density, group=1)) + # if only one data per date
  geom_line() +
  geom_smooth(se = TRUE) +
  labs(x = "Year", y = "Abundance") 


###################################################################################################################
Endgsdym <- read.csv("Abg_SDMdyEndog.csv") # for bar Animation
head(Endgsdym)
str(Endgsdym)
DFEndgDy <- data.frame(Endgsdym)
DFEndgDy <- na.omit(DFEndgDy)#Remove NA completely
DFEndgDy$X <- NULL# depend on the data tom remove unwanted column
head(DFEndgDy)
px <- ggplot(DFEndgDy, aes(Observation.date..Sampling.event.,  Density, fill =  Density)) +
  geom_col() +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "white"),
    panel.ontop = TRUE
  )
#animation
px + transition_states(Observation.date..Sampling.event., wrap = FALSE) +
  shadow_mark()
anim_save("Dyn_Endog.gif")

#data aggregation by date when multiple dates and then choose year interval e.g. 1970-1975
AggEndDy <- aggregate(DFEndgDy$Density, by=list(DFEndgDy$Observation.date..Sampling.event.), FUN=mean)
DF_AggEndDy <- data.frame(AggEndDy)
DF_AggEndDy <- na.omit(DF_AggEndDy)#Remove NA completely
head(DF_AggEndDy)

P <- ggplot(data = DF_AggEndDy, aes(x =Group.1, y = x, group=1)) + # if only one data per date
  geom_line() +
  geom_smooth(se = TRUE) +
  labs(x = "Year", y = "Abundance")




######################################################################################################################################################
library(dynamicSDM)#non working package
library(lubridate)# date conversion?

Endgsdym1 <- read.csv("DynamicSp_test.csv") # for bar Animation
head(Endgsdym1)
str(Endgsdym1)
DFEndgDy1 <- data.frame(Endgsdym1)
head(DFEndgDy1)
head(DFEndgDy1)
str(DFEndgDy1)
DFEndgDy1$year <- mdy(DFEndgDy1$year)
sample_occ_filtered <- spatiotemp_check(occ.data = DFEndgDy1 ,
                                        na.handle = "exclude",
                                        date.handle = "exclude",
                                        date.res = "year",
                                        coord.handle = "exclude",
                                        duplicate.handle ="exclude")
nrow(sample_occ_filtered)


bias_results <- spatiotemp_bias(occ.data = DFEndgDy1,
                                temporal.level = c("year"),
                                plot = FALSE,
                                spatial.method = "simple",
                                prj = "+proj=longlat +datum=WGS84")



DFEndgDy1$X <- NULL# depend on the data tom remove unwanted column
DFEndgDy1$X.1 <- NULL
DFEndgDy1$X.2 <- NULL
DFEndgDy1$X.3 <- NULL
DFEndgDy1$X.4 <- NULL
DFEndgDy1$X.5 <- NULL
DFEndgDy1$X.6 <- NULL
DFEndgDy1$X.7 <- NULL
DFEndgDy1$X.8 <- NULL
DFEndgDy1$X.9 <- NULL
DFEndgDy1$X.10 <- NULL
DFEndgDy1$X.11 <- NULL
DFEndgDy1$X.12 <- NULL
DFEndgDy1$X.13 <- NULL
DFEndgDy1$X.14 <- NULL
############################################################################################################################################################

EndgsdyCl <- read.csv("EndgCl_Dym2.csv") # for bar Animation
DF_endgCl_Dy <- data.frame (EndgsdyCl)
head(DF_endgCl_Dy)
#Aggregation
AggdyEndg <- aggregate(DF_endgCl_Dy$Ave..Abundance..site., by=list(DF_endgCl_Dy$Date, DF_endgCl_Dy$Habitat), FUN=mean)
AggdyEndg <- na.omit(AggdyEndg)#Remove NA completely
head(AggdyEndg, 20)
summary(AggdyEndg)
str(AggdyEndg)


px2 <- ggplot(AggdyEndg, aes(Group.1,   x, fill = x)) +
  geom_col() +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "white"),
    panel.ontop = TRUE
  )

px2 + transition_states(Group.1, wrap = TRUE) +
  shadow_mark()
#wrap
px2 +facet_wrap(~Group.2)+ transition_states(Group.1, wrap = TRUE) +
  shadow_mark()
anim_save("Dyn_Endog_wrap.gif")


################################################################sorted data
EndgsdySort <- read.csv("Sort_Dy Endog.csv")  # for bar graph Animation
DF_endgst_Dy <- data.frame (EndgsdySort)
head(DF_endgst_Dy)

px3 <- ggplot(DF_endgst_Dy, aes(Date,   Density , fill = Density )) +
  geom_col() +
  scale_fill_distiller(palette = "reds", direction = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "white"),
    panel.ontop = TRUE
  )

px3+facet_wrap(~Habitat)
#animated
px3 + transition_states(Date, wrap = TRUE) +
  shadow_mark()
anim_save("Dyn_Endog_sort.gif")# add gif to title

px3+facet_wrap(~Habitat) + transition_states(Date, wrap = TRUE) +
  shadow_mark()

anim_save("Dyn_Endog sort_wrap.gif")

#Line graph trends
DF_endgst_Dy[,'Date'] = as.factor(DF_endgst_Dy[,'Date'])#categorical to display indiv. date
Px4 <- ggplot(data = DF_endgst_Dy, aes(x =Date, y = Density , group=1)) + # if only one data per date
  geom_line() +
  geom_smooth(se = TRUE) +
  labs(x = "Year", y = "Abundance") 

Px4+ transition_reveal(Density)# animation not really meaningful

Px4+geom_point(aes(group = seq_along(Density))) +
  transition_reveal(Density)



