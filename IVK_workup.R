library(ggplot2)
library(reshape2)
library(xlsx)
library(plyr)
library(drc)
library(drfit)

#--- Data import and conversion to numeric ----
template <- read.xlsx("template2.xlsx", header = TRUE, sheetIndex = 1)
template <- template[-1,]

template$High.Concentration <- as.numeric(levels(template$High.Concentration)[template$High.Concentration])
template$High.Activity <- as.numeric(levels(template$High.Activity)[template$High.Activity])


#--- Set the number of ID columns ----
idCols <- 12


#--- Extrapolate concentration and activity data from initial values and serial dilution factor ----
conc <- matrix(nrow = nrow(template), ncol = max(template$num.wells))
for(i in seq(1:nrow(conc))){
  conc[i,1] = template$High.Concentration[i]
  for(j in seq(1:(template$num.wells[i]-1))+1){
    conc[i,j] <- conc[i,j-1] / template$Serial.Factor[i]
  }
}

act <- matrix(nrow = nrow(template), ncol = max(template$num.wells))
#for(i in seq(1:nrow(act))){
#act <- matrix(nrow = 8, ncol = 8)
for(i in seq(1:nrow(act))){
  act[i,1] = template$High.Activity[i]
  for(j in seq(1:(template$num.wells[i]-1))+1){
    act[i,j] <- act[i,j-1] / template$Serial.Factor[i]
  }
}


#--- Melt the data to pairs, add concentration and activity, cleanup ----
mtemp <- melt(template[,1:(idCols + nrow(conc))], id = colnames(template)[1:idCols])
mtemp <- mtemp[with(mtemp, order(ExpNum)), ]
mtemp <- cbind(mtemp, melt(t(conc))[,3])
mtemp <- cbind(mtemp, melt(t(act))[,3])
mtemp <- mtemp[, -(idCols+1)]


#--- Get the STD data from the template import, rename cols ----
mtemp <- cbind(mtemp, melt(t(template[,21:28]))[,3])
colnames(mtemp) <- c(colnames(mtemp)[1:12], 'Viability', 'Concentration', 'Activity', 'STD')


#---- Convert data and errors to percentages ----
mtemp$Viability <- mtemp$Viability * 100
mtemp$STD <- mtemp$STD * 100

#---- CURVE FITTING  ----
fits <- data.frame(ExpNum = 1, Concentration = 1, y = 1)
params <- matrix(nrow = nrow(template), ncol = 4)
drfun <- function(x, b, c, d, e) c + (d - c) / (1 + exp(b*(log(x) - log(e))))
for(i in seq(nrow(template))){
  rr <- drm(Viability~Concentration, 
            data = mtemp[ which(mtemp$ExpNum == i), ], 
            fct = LL.4())
  params[i, ] <- c(rr$fit$par[1], rr$fit$par[2], rr$fit$par[3], rr$fit$par[4])
  vals <- with(mtemp, seq(min(mtemp$Concentration[ which(mtemp$ExpNum == i)]), 
                          max(mtemp$Concentration[ which(mtemp$ExpNum == i)]), 
                          length = 500))
  addfits <- ddply(mtemp[ which(mtemp$ExpNum == i), ], "ExpNum", function(drfun) {
    data.frame(
      Concentration = vals, 
      y = drfun(vals, params[i, 1], 
                params[i, 2], 
                params[i, 3], 
                params[i, 4])
    )
  })
  fits <- rbind(fits, addfits)
}  
fits <- fits[-1, ]


#---------- END OF DATA MANIPULATION !-!-! PLOTTING STARTS ----------------
w = 0.65
fwid = 9
fhei = 6
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))


ggplot(mtemp, aes(x = Concentration, y = Viability, by = ExpNum)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), color="black", width=.1) +
  geom_point() +
  geom_line(aes(y = y), data = fits, color = "red") +
  scale_x_log10() +
  facet_wrap(~ ExpNum) +
  labs(x = "Concentration (nM)", y = "% Viability")


ggplot(mtemp, aes(x = Activity, y = Viability, by = ExpNum)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), color="black", width=.1) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~ExpNum) +
  labs(x = "Concentration (nM)", y = "% Viability")
