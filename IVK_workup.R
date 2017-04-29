library(ggplot2)
library(reshape2)
library(xlsx)
library(plyr)
library(drc)
library(drfit)

#--- Data import and conversion to level and numeric ----
template <- read.xlsx("template2.xlsx", header = TRUE, sheetIndex = 1)
template <- template[-1,]

template$High.Concentration <- as.numeric(levels(template$High.Concentration)[template$High.Concentration])
template$High.Activity <- as.numeric(levels(template$High.Activity)[template$High.Activity])
template$ExpNum <- as.factor(template$ExpNum)

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

#---- CURVE FITTING  SETUP ----
curvelength = 500
fits <- data.frame(ExpNum = 1, Concentration = 1, Activity = 1, Fit.Conc = 1, Fit.Act = 1)
cparams <- matrix(nrow = nrow(template), ncol = 4)
aparams <- matrix(nrow = nrow(template), ncol = 4)
drfun <- function(x, b, c, d, e) c + (d - c) / (1 + exp(b*(log(x) - log(e))))


for(i in seq(nrow(template))){
  #---- Build the dose response for Concentration ----
  crr <- drm(Viability~Concentration, 
            data = mtemp[ which(mtemp$ExpNum == i), ], 
            fct = LL.4())
  cparams[i, ] <- c(crr$fit$par[1], crr$fit$par[2], crr$fit$par[3], crr$fit$par[4])
  
  #---- Make the abscissa vector using curvelength ----
  vals <- with(mtemp, seq(min(mtemp$Concentration[ which(mtemp$ExpNum == i)]), 
                          max(mtemp$Concentration[ which(mtemp$ExpNum == i)]), 
                          length = curvelength))
  
  #---- Subset the mtemp frame by ExpNum and make the curve from fit params ----
  addfits <- ddply(mtemp[ which(mtemp$ExpNum == i), ], "ExpNum", function(drfun) {
    data.frame(
      Concentration = vals, 
      Fit.Conc = drfun(vals, cparams[i, 1], 
                cparams[i, 2], 
                cparams[i, 3], 
                cparams[i, 4])
    )
  })
  
  #---- Build the dose response for Activity ----
  arr <- drm(Viability~Activity, 
             data = mtemp[ which(mtemp$ExpNum == i), ], 
             fct = LL.4())
  
  #---- Make the abscissa vector using curvelength ----
  aparams[i, ] <- c(arr$fit$par[1], arr$fit$par[2], arr$fit$par[3], arr$fit$par[4])
  avals <- with(mtemp, seq(min(mtemp$Activity[ which(mtemp$ExpNum == i)]), 
                          max(mtemp$Activity[ which(mtemp$ExpNum == i)]), 
                          length = curvelength))
  
  #---- Subset the mtemp frame by ExpNum and make the curve from fit params ----
  aaddfits <- ddply(mtemp[ which(mtemp$ExpNum == i), ], "ExpNum", function(drfun) {
    data.frame(
      Activity = avals, 
      Fit.Act = drfun(avals, aparams[i, 1], 
                       aparams[i, 2], 
                       aparams[i, 3], 
                       aparams[i, 4])
    )
  })
  #---- Merge the concentration and activity fits ----
  addfits <- cbind(addfits, aaddfits[, -1])
  
  #---- Add each ExpNum to the end of fits
  fits <- rbind(fits, addfits)
}  

#---- Clean up the fits, add the categorical variables back in ----
fits <- fits[-1, ]
fitvars <- template[, 1:idCols]

#---- Expand, melt, and append the variables to the fits ----
fitvars <- cbind(fitvars, matrix(nrow = nrow(template), ncol = curvelength))
mfitvars <- melt(fitvars, id = colnames(template)[1:idCols])
fits <- cbind(fits[with(fits, order(ExpNum)), ], mfitvars[with(mfitvars, order(ExpNum)), 2:idCols])


#---------- END OF DATA MANIPULATION !-!-! PLOTTING STARTS ----------------
w = 0.65
fwid = 9
fhei = 6
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))


ggplot(mtemp, aes(x = Concentration, 
                  y = Viability, 
                  group = ExpNum, 
                  color = Cell.Line)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), width=.1) +
  geom_point() +
  geom_line(aes(y = Fit.Conc, linetype = Construct), data = fits) +
  facet_wrap(~ Antibody) +
  scale_x_log10() +
  labs(x = "Concentration (nM)", y = "% Viability")


ggplot(mtemp, aes(x = Activity, 
                  y = Viability, 
                  group = ExpNum, 
                  color = Cell.Line)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), width=.1) +
  geom_point() +
  geom_line(aes(y = Fit.Act, linetype = Construct), data = fits) +
  scale_x_log10() +
  facet_wrap(~ Antibody) +
  labs(x = "Activity (nCi)", y = "% Viability")
