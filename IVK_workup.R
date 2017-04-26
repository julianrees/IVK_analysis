library(ggplot2)
library(reshape2)
library(xlsx)
library(plyr)
rm(template)

template <- read.xlsx("template2.xlsx", header = TRUE, sheetIndex = 1)
template <- template[-1,]

template$High.Concentration <- as.numeric(levels(template$High.Concentration)[template$High.Concentration])
template$High.Activity <- as.numeric(levels(template$High.Activity)[template$High.Activity])


conc <- matrix(nrow = 8, ncol = 8)
for(i in seq(1:nrow(template))){
  conc[i,1] = template$High.Concentration[i]
  for(j in seq(1:7)+1){
    print(j)
    conc[i,j] <- conc[i,j-1] / template$Serial.Factor[i]
  }
}

act <- matrix(nrow = 8, ncol = 8)
for(i in seq(1:nrow(template))){
  act[i,1] = template$High.Activity[i]
  for(j in seq(1:7)+1){
    print(j)
    act[i,j] <- act[i,j-1] / template$Serial.Factor[i]
  }
}
melt(t(conc))[,3]
rm(mtemp)
mtemp <- melt(template[,1:20], id = colnames(template)[1:idCols])
mtemp <- mtemp[with(mtemp, order(ExpNum)), ]
mtemp <- cbind(mtemp, melt(t(conc))[,3])
mtemp <- cbind(mtemp, melt(t(act))[,3])

# - - - STD


mtemp <- cbind(mtemp, melt(t(template[,21:28]))[,3])
drop(mtemp$variable)


colnames(mtemp) <- c(colnames(mtemp)[1:12], 'n1', 'Viability', 'Concentration', 'Activity', 'STD')

mtemp$Viability <- mtemp$Viability * 100
mtemp$STD <- mtemp$STD * 100

#---------- END OF DATA MANIPULATION !-!-! PLOTTING STARTS ----------------
w = 0.65
fwid = 9
fhei = 6
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))


ggplot(mtemp, aes(x = Concentration, y = Viability, by = ExpNum)) + 
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), color="black", width=.1) + 
  geom_point() + 
  scale_x_log10() + 
  facet_wrap(~ExpNum) + 
  labs(x = "Concentration (nM)", y = "% Viability")
  

ggplot(mtemp, aes(x = Activity, y = Viability, by = ExpNum)) + 
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), color="black", width=.1) + 
  geom_point() + 
  scale_x_log10() + 
  facet_wrap(~ExpNum)
  












