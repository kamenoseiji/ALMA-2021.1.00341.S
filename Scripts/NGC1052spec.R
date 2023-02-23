#setwd('/Users/kameno/Dropbox/GBT_NGC1052')
#setwd('/Volumes/SSD/GBT_NGC1052')
FWHM_sigma <- 2.0* sqrt(2.0*log(2.0))
textVel = 980  # Label position (km/s)
fillColor <- c(rgb(1, 0.75, 0.75), rgb(0.75, 1, 0.75), rgb(0.75, 0.75, 1))
#setwd('/home/skameno/GBT_NGC1052')
source('GBTfunction.R')
lumiDistance <- 17.6    # Kameno+2020
luminosity <- function(Freq, Flux, Distance){                     # Hagiwara+2021 ApJ 923, 251
    # Freq : rest frequency [GHz]
    # Flux : integrated flux density [Jy km/s]
    # Distance : luminosity distance [Mpc]
    return(1.04e-3 * Freq* Flux* Distance^2)     # in unit of solar luminosity
}
fillSpec <- function(veloc, flux, baseLevel, color){       # fill the spectrum with color
    chWid <- median(diff(veloc))
    plotVeloc <- c(veloc[1]- 0.5*chWid,  as.vector(t(matrix(c(veloc - 0.5*chWid, veloc + 0.5*chWid), ncol=2))), veloc[length(veloc)] + 0.5*chWid)
    plotFlux  <- c(baseLevel, as.vector(t(matrix(rep(flux, 2), ncol=2))) + baseLevel, baseLevel)
    polygon(plotVeloc, plotFlux, col=color, border=NA)
}

pdf("NGC1052_maser_series.pdf", width=7, height=11)


#-------- 22-GHz plots
fileList <- c('AGBT05C_034_02', 'AGBT05C_034_03', 'AGBT05C_034_04', 'AGBT05C_034_05', 'AGBT05C_034_07', 'AGBT05C_034_10')
dateList <- c('2005-11-18', '2005-11-28', '2005-12-07', '2005-12-20', '2005-12-27', '2006-01-13')
lineLabel <- expression('22 GHz H'[2]*'O 6'[1][{","}][6]* '-5'[2][{","}][3])
weightList <- list(
    c(rep(0, 263), rep(1,(2563-263)), rep(0,(5100-2563)), rep(1,(7975-5100)), rep(0,(8192-7975))),   # AGBT05C_034_02
    c(rep(0, 263), rep(1,(2563-263)), rep(0,(4585-2563)), rep(1,(7975-4585)), rep(0,(8192-7975))),   # AGBT05C_034_03
    c(rep(0, 263), rep(1,(2563-263)), rep(0,(4585-2563)), rep(1,(7975-4585)), rep(0,(8192-7975))),   # AGBT05C_034_04
    c(rep(0, 263), rep(1,(2563-263)), rep(0,(4585-2563)), rep(1,(7975-4585)), rep(0,(8192-7975))),   # AGBT05C_034_05
    c(rep(0, 239), rep(1,(1156-239)), rep(0,(2709-1156)), rep(1,(7964-2709)), rep(0,(8192-7964))),   # AGBT05C_034_07
    c(rep(0, 212), rep(1,(829-212)), rep(0,(2597-829)), rep(1,(7939-2597)), rep(0,(8192-7939))))     # AGBT05C_034_10
plot(weightList[[1]], type='n', axes=TRUE, xlim=c(1000,2000), ylim=c(-0.00,0.6), xlab='LSR Velocity [km/s]', ylab='Flux Density [Jy]')

fillRange <- c(1200, 2000)
for(file_index in 1:6){
    levelBase <- 0.05*(10.5-file_index)
    df <- fitdata(baselinefit(sprintf('%s.txt', fileList[file_index]), weightList[[file_index]]))
    chWid <- median(diff(df$velocity))
    columnPlotVel <- df$velocity - 0.5*chWid
    fillIndex <- which( (df$velocity > fillRange[1]) & (df$velocity < fillRange[2]))
    fillSpec(df$velocity[fillIndex], df$flux[fillIndex], levelBase, fillColor[1])
    lines(df$velocity - 0.5*chWid, df$flux + levelBase, type='s')
    abline(h=levelBase, col='gray60', lwd=0.1)
    #---- calculate integrated flux density and luminosity
    DF <- df[((df$velocity > 1400) & (df$velocity < 1900)),]    # trim the velocity range
    integFlux <- sum(DF$flux) * abs(median(diff(DF$velocity)))  # [Jy km/s]
    isoLumi <- luminosity(22.235, integFlux, lumiDistance)
    text_sd <- sprintf('Flux: %.1f Jy km/s  Luminosity: %.1f solar luminosity\n', integFlux, isoLumi)
    cat(text_sd)
    text(textVel, levelBase + 0.01, lineLabel, cex=0.7, adj=0, pos=4, col='darkred')
    text(textVel, levelBase + 0.02, dateList[file_index], cex=0.7, adj=0, pos=4, col='darkred')
}

#---- HCN (1-0) absorption
levelBase <- 0.2
lineLabel <- '89 GHz HCN 1-0'
fillRange <- c(1120, 1800)
df <- read.table('HCN_J=1-0.data', header=T)
chWid <- median(diff(df$VLSR))
columnPlotVel <- df$VLSR - 0.5*chWid
fillIndex <- which( (df$VLSR > fillRange[1]) & (df$VLSR < fillRange[2]))
fillSpec(df$VLSR[fillIndex], (df$flux[fillIndex] - df$refSpec[fillIndex]), levelBase, fillColor[2])
lines(df$VLSR - 0.5*chWid, (df$flux - df$refSpec) + levelBase, type='s')
abline(h=levelBase, col='gray60', lwd=0.1)
text_sd = '2017-08-28'
text(textVel, levelBase+0.005, paste(lineLabel, text_sd, sep=' '), cex=0.7, adj=0, pos=4, col='darkgreen')
# text(textVel, levelBase+0.005, lineLabel, cex=0.7, adj=0, pos=4, col='darkgreen')
# text(textVel, levelBase+0.02, text_sd, cex=0.7, adj=0, pos=4, col='darkgreen')

#---- H2O_5(5,0)-6(4,3)_v2=1 absorption
levelBase <- 0.16
lineLabel <- expression('233 GHz H'[2]*'O 5'[5][{","}][0]* '-6'[4][{","}][3]* ' x15')
fillRange <- c(1200, 1900)
df <- read.table('H2O_5(5,0)-6(4,3)_v2=1.data', header=T)
chWid <- median(diff(df$VLSR))
columnPlotVel <- df$VLSR - 0.5*chWid
fillIndex <- which( (df$VLSR > fillRange[1]) & (df$VLSR < fillRange[2]))
fillSpec(df$VLSR[fillIndex], 15*(df$flux[fillIndex] - df$refSpec[fillIndex]), levelBase, fillColor[2])
lines(df$VLSR - 0.5*chWid, 15*(df$flux - df$refSpec) + levelBase, type='s')
abline(h=levelBase, col='gray60', lwd=0.1)
text_sd = '2015-08-05'
text(textVel, levelBase + 0.005, lineLabel, cex=0.7, adj=0, pos=4, col='darkgreen')
text(textVel, levelBase + 0.015, text_sd, cex=0.7, adj=0, pos=4, col='darkgreen')


#par(mfrow = c(2, 1))
#par(oma=c(0,0,0,0), mar=(0,0,0,0))
#---- sub-mm maser
levelBase <- 0.0
lineLabel <- expression('321 GHz H'[2]*'O 10'[2][{","}][9]* '-9'[3][{","}][6])
fillRange <- c(1200, 1800)
fileList <- c('NGC1052_SPW2', 'NGC1052_SPW3')
df <- read.table(sprintf('%s.veloc.txt', fileList[1]), skip=8)
df <- rbind(df, read.table(sprintf('%s.veloc.txt', fileList[2]), skip=8))
names(df) <- c("velocity", "flux")
df$flux <- df$flux - median(df$flux)
df <- df[order(df$velocity),]
chWid <- c(diff(df$velocity), median(diff(df$velocity)))
#---- sub-mm maser decomposition
fit <- nls( formula = flux ~ a* exp(-0.5*((velocity - b)/c)^2) + d* exp(-0.5*((velocity - e)/f)^2), data=df, start=list(a=0.036, b=1468, c=88, d=0.05, e=1330, f=20))
text_sd <- sprintf('Comp1 : %.1f +- %.1f mJy : VLSR=%.1f +- %.1f  FWHM=%.1f +- %.1f\n',
                summary(fit)$coefficients[1]*1e3, summary(fit)$coefficients[7]*1e3,
                summary(fit)$coefficients[2], summary(fit)$coefficients[8],
                summary(fit)$coefficients[3]*FWHM_sigma, summary(fit)$coefficients[9]*FWHM_sigma)
cat(text_sd)
text_sd <- sprintf('Comp2 : %.1f +- %.1f mJy : VLSR=%.1f +- %.1f  FWHM=%.1f +- %.1f\n',
                summary(fit)$coefficients[4]*1e3, summary(fit)$coefficients[10]*1e3,
                summary(fit)$coefficients[5], summary(fit)$coefficients[11],
                summary(fit)$coefficients[6]*FWHM_sigma, summary(fit)$coefficients[12]*FWHM_sigma)
cat(text_sd)
fillIndex <- which( (df$velocity > fillRange[1]) & (df$velocity < fillRange[2]))
fillSpec(df$velocity[fillIndex], df$flux[fillIndex], levelBase, fillColor[3])
lines( df$velocity, fitted(fit) + levelBase, col='blue', lwd=2)
lines( df$velocity, summary(fit)$coefficients[1]* exp(-0.5*((df$velocity - summary(fit)$coefficients[2])/summary(fit)$coefficients[3])^2) + levelBase, col='darkgreen', lty=2, lwd=1.5)
lines( df$velocity, summary(fit)$coefficients[4]* exp(-0.5*((df$velocity - summary(fit)$coefficients[5])/summary(fit)$coefficients[6])^2) + levelBase, col='deeppink', lty=2, lwd=1.5)
#---- plot sub-mm maser actual data
lines(df$velocity - 0.5*chWid, df$flux + levelBase, type='s')
abline(h=levelBase, col='gray60', lwd=0.1)
#---- calculate integrated flux density and luminosity
DF <- df[((df$velocity > 1250) & (df$velocity < 1700)),]    # trim the velocity range
integFlux <- sum(DF$flux) * abs(median(diff(DF$velocity)))  # [Jy km/s]
restFreq <- 321.226
isoLumi <- luminosity(restFreq, integFlux, lumiDistance)
text_sd <- sprintf('Flux: %.1f Jy km/s  Luminosity: %.1f solar luminosity\n', integFlux, isoLumi)
cat(text_sd)
text_sd <- '2022-05-12'
text(textVel, levelBase + 0.02, lineLabel, cex=0.7, adj=0, pos=4, col='navy')
text(textVel, levelBase + 0.03, text_sd, cex=0.7, adj=0, pos=4, col='navy')


#---- HCN (4-3) absorption
#par(oma=c(0,0,0,0), mar=(0,0,0,0))
#plot(weightList[[1]], type='n', axes=TRUE, xlim=c(1000,2000), ylim=c(-0.06,0.0), xlab='LSR Velocity [km/s]', ylab='')
levelBase <- 0.11
lineLabel <- '355 GHz HCN 4-3'
fillRange <- c(1150, 1900)
# df <- read.table('HCN_J=4-3_v=0.data', header=T)
df <- read.table('HCN43_veloc.txt', header=F, skip=7)
names(df) <- c('VLSR', 'flux')
df$VLSR <- rev(df$VLSR)
df$flux <- rev(df$flux) - quantile(df$flux, 0.8)
chWid <- median(diff(df$VLSR))
columnPlotVel <- df$VLSR - 0.5*chWid
fillIndex <- which( (df$VLSR > fillRange[1]) & (df$VLSR < fillRange[2]))
fillSpec(df$VLSR[fillIndex], df$flux[fillIndex], levelBase, fillColor[2])
lines(df$VLSR - 0.5*chWid, df$flux + levelBase, type='s')
abline(h=levelBase, col='gray60', lwd=0.1)
text_sd = '2015-08-16'
text(textVel, levelBase+0.005, paste(lineLabel, text_sd, sep=' '), cex=0.7, adj=0, pos=4, col='darkgreen')
#text(textVel, levelBase+0.02, text_sd, cex=0.7, adj=0, pos=4, col='darkgreen')
abline(v=1492, col='gray60', lwd=0.2, lty=2)
dev.off()

Tb <- function(flux, freq, FWHM){
    # flux : [Jy]
    # freq : [GHz]
    # FWHM : [arcsec]
    return( 1.222e6* flux / freq^2 / FWHM^2 )
}