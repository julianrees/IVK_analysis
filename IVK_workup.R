library(ggplot2)
library(reshape2)
library(xlsx)
library(plyr)
library(drc)
library(drfit)

#--- Data import and conversion to level and numeric ----
experiments <- read.xlsx("template2.xlsx", header = TRUE, sheetIndex = 1)
experiments <- experiments[-1,]

experiments$High.Concentration <- as.numeric(levels(experiments$High.Concentration)[experiments$High.Concentration])
experiments$High.Activity <- as.numeric(levels(experiments$High.Activity)[experiments$High.Activity])
experiments$ExpNum <- as.factor(experiments$ExpNum)

#--- Set the number of ID columns ----
idCols <- 12


#--- Extrapolate concentration and activity data from initial values and serial dilution factor ----
conc <- matrix(nrow = nrow(experiments), ncol = max(experiments$num.wells))
for(i in seq(1:nrow(conc))){
  conc[i,1] = experiments$High.Concentration[i]
  for(j in seq(1:(experiments$num.wells[i]-1))+1){
    conc[i,j] <- conc[i,j-1] / experiments$Serial.Factor[i]
  }
}

act <- matrix(nrow = nrow(experiments), ncol = max(experiments$num.wells))
for(i in seq(1:nrow(act))){
  act[i,1] = experiments$High.Activity[i]
  for(j in seq(1:(experiments$num.wells[i]-1))+1){
    act[i,j] <- act[i,j-1] / experiments$Serial.Factor[i]
  }
}


#--- Melt the data to pairs, add concentration and activity, cleanup ----
mexp <- melt(experiments[,1:(idCols + nrow(conc))], id = colnames(experiments)[1:idCols])
mexp <- mexp[with(mexp, order(ExpNum)), ]
mexp <- cbind(mexp, melt(t(conc))[,3])
mexp <- cbind(mexp, melt(t(act))[,3])
mexp <- mexp[, -(idCols+1)]

#test cahnge

#--- Get the STD data from the experiments import, rename cols ----
mexp <- cbind(mexp, melt(t(experiments[,21:28]))[,3])
colnames(mexp) <- c(colnames(mexp)[1:12], 'Viability', 'Concentration', 'Activity', 'STD')


#---- Convert data and errors to percentages ----
#mexp$Viability <- mexp$Viability * 100
#mexp$STD <- mexp$STD * 100

#---- CURVE FITTING  SETUP ----
curvelength = 5000
fits <- data.frame(ExpNum = 1, Concentration = 1, Activity = 1, Fit.Conc = 1, Fit.Act = 1)
cparams <- matrix(nrow = nrow(experiments), ncol = 4)
aparams <- matrix(nrow = nrow(experiments), ncol = 4)
results <- data.frame(nrow = nrow(experiments), ncol = 4)
drfun <- function(x, b, c, d, e) c + (d - c) / (1 + exp(b*(log(x) - log(e))))

i = 1
for(i in seq(nrow(experiments))){
  #---- Build the dose response for Concentration ----
  crr <- drm(Viability~Concentration,
            data = mexp[ which(mexp$ExpNum == i), ],
            fct = LL.4())
  cparams[i, ] <- c(crr$fit$par[1], crr$fit$par[2], crr$fit$par[3], crr$fit$par[4])
  
  #---- Make and fit the low and high error curves
  low <- mexp[ which(mexp$ExpNum == i), ]
  high <- mexp[ which(mexp$ExpNum == i), ]
  low$Viability <- low$Viability - low$STD
  high$Viability <- high$Viability + low$STD
  
  lccr <- drm(Viability~Concentration, data = low, fct = LL.4())
  hccr <- drm(Viability~Concentration, data = high, fct = LL.4())
  
  #---- Make the abscissa vector using curvelength ----
  vals <- with(mexp, seq(min(mexp$Concentration[ which(mexp$ExpNum == i)]),
                          max(mexp$Concentration[ which(mexp$ExpNum == i)]),
                          length = curvelength))

  #---- Subset the mexp frame by ExpNum and make the curve from fit params ----
  addfits <- ddply(mexp[ which(mexp$ExpNum == i), ], "ExpNum", function(drfun) {
    data.frame(
      Concentration = vals,
      Fit.Conc = drfun(vals, cparams[i, 1],
                cparams[i, 2],
                cparams[i, 3],
                cparams[i, 4])
    )
  })

  #---- Put the results in to table
  results[i,1] <- crr$fit$par[4]*1000
  results[i,2] <- lccr$fit$par[4]*1000
  results[i,3] <- hccr$fit$par[4]*1000

  #---- Build the dose response for Activity ----
  arr <- drm(Viability~Activity,
             data = mexp[ which(mexp$ExpNum == i), ],
             fct = LL.4())

  #---- Make and fit the low and high error curves
  larr <- drm(Viability~Activity, data = low, fct = LL.4())
  harr <- drm(Viability~Activity, data = high, fct = LL.4())
  
  
  #---- Make the abscissa vector using curvelength ----
  aparams[i, ] <- c(arr$fit$par[1], arr$fit$par[2], arr$fit$par[3], arr$fit$par[4])
  avals <- with(mexp, seq(min(mexp$Activity[ which(mexp$ExpNum == i)]),
                          max(mexp$Activity[ which(mexp$ExpNum == i)]),
                          length = curvelength))

  #---- Subset the mexp frame by ExpNum and make the curve from fit params ----
  aaddfits <- ddply(mexp[ which(mexp$ExpNum == i), ], "ExpNum", function(drfun) {
    data.frame(
      Activity = avals,
      Fit.Act = drfun(avals, aparams[i, 1],
                       aparams[i, 2],
                       aparams[i, 3],
                       aparams[i, 4])
    )
  })

  #---- Put the results in to table
  results[i,4] <- arr$fit$par[4]*1000
  results[i,5] <- larr$fit$par[4]*1000
  results[i,6] <- harr$fit$par[4]*1000
  
  #---- Merge the concentration and activity fits ----
  addfits <- cbind(addfits, aaddfits[, -1])

  #---- Add each ExpNum to the end of fits
  fits <- rbind(fits, addfits)
}

#---- Clean up the fits, add the categorical variables back in ----
fits <- fits[-1, ]
fitvars <- experiments[, 1:idCols]

#---- Expand, melt, and append the variables to the fits ----
fitvars <- cbind(fitvars, matrix(nrow = nrow(experiments), ncol = curvelength))
mfitvars <- melt(fitvars, id = colnames(experiments)[1:idCols])
fits <- cbind(fits[with(fits, order(ExpNum)), ], mfitvars[with(mfitvars, order(ExpNum)), 2:idCols])

#---- Add results into the experiment table
colnames(results) <- c('Concentration.EC50', 'Concentration.LowEC50Limit','Concentration.HighEC50Limit', 
                       'Activity.EC50', 'Activity.LowEC50Limit', 'Activity.HighEC50Limit')
experiments <- cbind(experiments, results)


#---------- END OF DATA MANIPULATION !-!-! PLOTTING STARTS ----------------
w = 0.65
fwid = 9
fhei = 6
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())


#mexp <- subset(mexp, Targeted == TRUE)
#fits <- subset(fits, Targeted == TRUE)


ggplot(mexp, aes(x = Concentration,
                  y = Viability,
                  group = ExpNum,
                  color = Antibody)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), width=.1) +
  geom_point(aes(shape = Construct), size = 3) +
  geom_line(aes(y = Fit.Conc, linetype = Construct), data = fits) +
  facet_wrap(~ ExpNum) +
  scale_x_log10(limits =c(min(mexp$Concentration),max(mexp$concentration))) +
  scale_y_continuous(breaks = seq(0,1,length.out=11), labels = scales::percent) +
  labs(x = "Concentration (nM)", y = "Viability")



ggplot(mexp, aes(x = Activity,
                  y = Viability,
                  group = ExpNum,
                  color = Cell.Line)) +
  geom_errorbar(aes(ymin=Viability-STD, ymax=Viability+STD), width=.1) +
  geom_point(aes(shape = Construct), size = 3) +
  geom_line(aes(y = Fit.Act, linetype = Construct), data = fits) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0,1,length.out=11), labels = scales::percent) +
  facet_wrap(~ Antibody) +
  labs(x = "Activity (nCi)", y = "% Viability")
  #change the variables above!
  ggsave(filename = 'Rfigs/thisthang.png',
       width = fwid, height = fhei, units = "in")


ggplot(experiments, aes(x = ExpNum, y = Concentration.EC50)) +
  geom_col(aes(fill = Antibody, color = Antibody, alpha = Cell.Line)) +
  geom_errorbar(aes(ymin = Concentration.LowEC50Limit,
                    ymax = Concentration.HighEC50Limit)) +
  scale_y_log10() +
  labs(x = "Experiment Number", y = expression(paste("IC"[50], " (pM)")))

