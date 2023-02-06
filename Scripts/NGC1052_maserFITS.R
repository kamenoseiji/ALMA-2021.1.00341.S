library(FITSio)
setwd('~/Projects/PRTSPR-58044')
pdf('NGC1052_Cont+Maser_Map.pdf', width=8, height=8)
pixSize <- 0.1  # 0.1 arcsec/pixel
ContFile <- 'NGC1052.ContUSB.fits'
LineFile <- 'NGC1052.LineCH.fits'

LineFITS <- readFITS(LineFile)
ContFITS <- readFITS(ContFile)

numXpix <- dim(LineFITS$imDat)[1]   # X-axis pixel
numYpix <- dim(LineFITS$imDat)[2]   # X-axis pixel

chColor <- rev(rainbow(21)[1:16])
VLSR <- c(
1280,
1300,
1320,
1340,
1360,
1380,
1400,
1420,
1440,
1460,
1480,
1500,
1520,
1540,
1560,
1580)

velLabel <- character(length=16)
contMap <- ContFITS$imDat[(ceiling(numXpix/2)+16):(ceiling(numXpix/2)-15), (ceiling(numYpix/2)-15):(ceiling(numYpix/2)+16) ,1,1]
RA <- DEC <- 0.1*(-16:15)
image(RA, DEC, contMap, col=grey.colors(256, start=-0.1, end=1.0, gamma=2), xlab='Relative RA [arcsec]', ylab='Relative Dec [arcsec]', xlim=c(1.5, -1.6), useRaster=TRUE)

numCH <- dim(LineFITS$imDat)[3]     # Number of channels
Xpos <- Ypos <- Xerr <- Yerr <- peakFlux <- numeric(numCH)
numXpix <- dim(LineFITS$imDat)[1]   # X-axis pixel
numYpix <- dim(LineFITS$imDat)[2]   # X-axis pixel
for(chIndex in 1:numCH){            # Loop for channel
    chMap <- data.frame(flux=as.vector(LineFITS$imDat[(ceiling(numXpix/2)-15):(ceiling(numXpix/2)+16), (ceiling(numYpix/2)-15):(ceiling(numYpix/2)+16) ,chIndex,1]), X=rep(16:(-15), 32), Y=as.vector(t(matrix(rep(-16:15, 32), nrow=32))))
    fit <- nls(data=chMap, formula=flux ~ a* exp(-0.5*((X-c)^2 + (Y-d)^2) / b^2 ), start=list(a=max(chMap$flux), b=3, c=0, d=0))
    Xpos[chIndex] <- pixSize*summary(fit)$coefficients[3]
    Ypos[chIndex] <- pixSize*summary(fit)$coefficients[4]
    Xerr[chIndex] <- 3.0*pixSize*summary(fit)$coefficients[7]
    Yerr[chIndex] <- 3.0*pixSize*summary(fit)$coefficients[8]
    peakFlux[chIndex] <- summary(fit)$coefficients[1]
    text_sd <- sprintf('CH%02d : X=%.2f +- %.3f Y=%.2f +- %.3f Peak=%.5f (%f) \n', chIndex, Xpos[chIndex], Xerr[chIndex], Ypos[chIndex], Yerr[chIndex], peakFlux[chIndex],  peakFlux[chIndex]/sd(chMap$flux))
    velLabel[chIndex] <- sprintf('%d km/s', VLSR[chIndex])
    cat(text_sd)
}

arrows(Xpos-Xerr, Ypos, Xpos+Xerr, Ypos, length=0, col=chColor)
arrows(Xpos, Ypos-Yerr, Xpos, Ypos+Yerr, length=0, col=chColor)
points(Xpos, Ypos, pch=20, col=chColor, cex=25*peakFlux)
legend("topright", legend=velLabel, col=chColor, pch=20, bg='white')

#-------- Scale bar
Scale10pc <- 0.1176471  # 0.1176471 arcsec corresponds to 10 pc
arrows(-0.5, -1.4, -0.5+Scale10pc, -1.4, length=0, col='white', lwd=2)
text(-0.45, -1.48, '10 pc', col='white')
dev.off()