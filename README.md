# ALMA-2021.1.00341.S
Data and processing script for ALMA 2021.1.00341.S (only for phase calibrator and bandpass scans)

Data reduction process and results from ALMA 2021.1.00341.S (only for phase calibrator and bandpass scans).

Scripts/ : reduction scripts (Python and R)
- NGC1052.py     : CASA data reduction script
- UVFIT.py       : CASA script to estimate the size using visibilities as a function of spatial frequencies (projected baselines)
- NGC1052spec.R  : R script to plot continuum-subracted spectra
- GBTfunction.R  : Velocity decomposition, called from NGC1052spec.R
- NGC1052_maserFITS.R : Plot continuum (contour) map and maser positions using FITS image files

Data/  : 
- J0241-0815_ContUSB.ms.xz : compressed MS file for continuum visibilities
- J0241-0815_Maser.ms.xz   : compressed MS file for continuum-subtracted visibilities
- NGC1052.ContUSB.fits     : FITS image file of continuum
- NGC1052.LineCH.fits      : FITS image file of continuum-subtracted line emission
- *.txt                    : Spectral data in ASCII-text format.

