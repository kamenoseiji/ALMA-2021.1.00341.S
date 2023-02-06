SCR_DIR = '/Users/skameno/ALMA_Py3/'
exec(open(SCR_DIR + 'interferometry.py').read())
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#-------- Gaussian visibility model 
def GaussVis(uvwave, theta):    # 1-dimensional Gaussian visibility (complex)
    # uvwave : uv distance in the unit of wavelength
    # theta  : FWHM in arcsec
    Factor = np.pi* (0.5 + 0.0j) / 206264.80624709636 / np.sqrt(np.log(2.0))
    return np.exp( -(Factor* theta* uvwave)**2 )
#
def GaussResid(vis, flux, uvwave, theta): # Calculate residual = observed visisbility - Gaussian model
    residVis = vis - flux* GaussVis(uvwave, theta)
    return  abs(residVis.dot(residVis.conjugate()))
#
def quadMin(x, y):
    a = y[0]/(x[0] - x[1])/(x[0] - x[2]) + y[1]/(x[1] - x[2])/(x[1] - x[0]) + y[2]/(x[2] - x[0])/(x[2] - x[1])
    b = y[0]* (x[1] + x[2])/(x[0] - x[1])/(x[0] - x[2]) + y[1]* (x[2] + x[0])/(x[1] - x[2])/(x[1] - x[0]) + y[2]* (x[0] + x[1])/(x[2] - x[0])/(x[2] - x[1])
    return 0.5*abs(b/a)
#
def visFit(spatialFreq, Vis):
    flux, size, step = abs(np.mean(Vis)), 0.1, 0.5
    residImag = 2.0* abs(Vis.imag.dot(Vis.imag))
    dof       = len(Vis) - 2
    while( step > 0.001):
        # print('flux=%.3f size=%e resid=%.1f/%d ' % (flux, size, dof*GaussResid(Vis, flux, spatialFreq, size)/residImag,dof))
        fluxTemp = np.array([(1.0+step)*flux, flux, flux/(1.0+step)])
        sizeTemp = np.array([(1.0+step)*size, size, size/(1.0+step)])
        resid = [GaussResid( Vis, flux, spatialFreq, sizeTemp[0]),
                 GaussResid( Vis, flux, spatialFreq, sizeTemp[1]),
                 GaussResid( Vis, flux, spatialFreq, sizeTemp[2])]
        if np.std(resid) < 0.001: break
        size = quadMin(sizeTemp,resid)
        resid = [GaussResid( Vis, fluxTemp[0], spatialFreq, size),
                 GaussResid( Vis, fluxTemp[1], spatialFreq, size),
                 GaussResid( Vis, fluxTemp[2], spatialFreq, size)]
        if np.std(resid) < 0.001: break
        flux = quadMin(fluxTemp,resid)
        step = 0.75* step
    #
    bestSize, bestFlux, bestResid = size, flux, dof*GaussResid(Vis, flux, spatialFreq, size)/residImag
    step = bestSize * 0.001
    minSize, maxSize = 0.0, bestSize
    minResid, maxResid = dof*GaussResid(Vis, bestFlux, spatialFreq, minSize)/residImag, bestResid
    #---- minimum acceptable size
    while minResid > bestResid + 2.57: # 99% confidence level
        minSize = minSize + step
        minResid = dof*GaussResid(Vis, bestFlux, spatialFreq, minSize)/residImag
    #---- maximum acceptable size
    while maxResid < bestResid + 2.57: # 99% confidence level
        maxSize = maxSize + step
        maxResid = dof*GaussResid(Vis, bestFlux, spatialFreq, maxSize)/residImag
    #
    return bestFlux, bestSize, minSize, maxSize
#
#-------- Maser and Cont data
maserMS = 'J0241-0815_Maser.ms'
contMS  = 'J0241-0815_ContUSB.ms'
scanList = [6,10,14,17]
chNum, chWid, freq = GetChNum(maserMS, 0)
maserVisList, contVisList, uvList = [], [], []
for scan in scanList:
    timeStamp, uvw = GetUVW(maserMS, 0, scan)
    timeStamp, maserPS, maserVis =  GetVisAllBL(maserMS, 0, scan)
    timeStamp, contPS,  contVis  =  GetVisAllBL(contMS, 0, scan)
    uvw = np.mean(uvw, axis=2)
    uvList = uvList + [np.sqrt(uvw[0]**2 + uvw[1]**2)]
    maserVisList = maserVisList + [np.mean(maserVis, axis=(0,1,3))]
    contVisList  = contVisList  + [np.mean(contVis, axis=(0,1,3))]
#
scanNum = len(scanList)
blNum = uvw.shape[1]
uvWave = np.array(uvList).reshape(scanNum*blNum)* freq[0] / 299792458   # uv distance in unit of wavelength
maserVis = np.array(maserVisList).reshape(scanNum*blNum)
contVis  = np.array(contVisList).reshape(scanNum*blNum)
flagIndex = np.where(abs(contVis) > 0.1* np.median(abs(contVis)))[0]
useBLnum = int(len(flagIndex)/scanNum)
plt.plot(np.mean(uvWave[flagIndex].reshape(scanNum, useBLnum), axis=0), abs(np.mean(contVis[flagIndex].reshape(scanNum, useBLnum), axis=0)), 'k.', label='continuum')
plt.plot(np.mean(uvWave[flagIndex].reshape(scanNum, useBLnum), axis=0), abs(np.mean(maserVis[flagIndex].reshape(scanNum, useBLnum), axis=0)), 'b.', label='H$_2$O line')
plt.xlabel('Spatial frequency [λ]')
plt.ylabel('Visibility amplitude [Jy]')
#-------- Continuum Size
contFlux, contSize, minContSize, maxContSize = visFit( uvWave[flagIndex], contVis[flagIndex] )
maserFlux, maserSize, minMaserSize, maxMaserSize = visFit( uvWave[flagIndex], maserVis[flagIndex] )
modelUVwave = np.arange(0.0, max(uvWave[flagIndex]), 100)
plt.fill_between(modelUVwave, contFlux* abs(GaussVis(modelUVwave, minContSize)), contFlux* abs(GaussVis(modelUVwave, maxContSize)), facecolor='black', alpha=0.25)
plt.fill_between(modelUVwave, maserFlux* abs(GaussVis(modelUVwave, minMaserSize)), maserFlux* abs(GaussVis(modelUVwave, maxMaserSize)), facecolor='navy', alpha=0.25)
print('Continuum: flux=%.3f bestSize=%e maxSize=%e' % (contFlux, contSize, maxContSize))
print('H$_2$O Line: flux=%.3f bestSize=%e maxSize=%e' % (maserFlux, maserSize, maxMaserSize))
plt.legend()
plt.savefig('UVplot.pdf')
plt.close('all')
