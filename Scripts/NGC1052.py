prefix =  'uid___A002_Xf8d822_X9d6a'
from almahelpers_localcopy import tsysspwmap
'''
#Fields: 3
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0423-0120          04:23:15.800720 -01.20.33.06550 ICRS    0       13837623
#  1    none J0241-0815          02:41:04.798502 -08.15.20.75189 ICRS    1        3154909
#  2    none XMM3-3085           02:26:59.080000 -05.12.17.49000 ICRS    2       18501865
'''
BPCAL = 'J0423-0120'
PHCAL= 'J0241-0815'
TARGET= 'XMM3-3085'
REFANT = 'DA52'
WVRSPW = '4'
'''
os.system('rm -rf ' + prefix + '_flagonline.txt')
importasdm(prefix)
plotants( vis=prefix + '.ms', figfile=prefix + '_plotants.png')
browsetable(tablename = prefix +'.ms')
listobs(vis=prefix+'.ms', scan='', spw='', verbose=True, listfile=prefix+'.listobs')
msfile = prefix + '.ms'
flagdata(vis=msfile,mode='manual', autocorr=True, flagbackup=False)
flagdata(vis=msfile,mode='manual', intent='*POINTING*,*ATMOSPHERE*', flagbackup=False)
flagdata(vis=msfile,mode='manual', antenna='DV05', flagbackup=False)
flagmanager(vis=msfile,mode='save',versionname='Original')
#---- Tsys table
os.system('rm -rf *T0.tsys')
gencal(vis = msfile, caltable = prefix + 'T0.tsys', caltype = 'tsys')
plotms(prefix + 'T0.tsys', xaxis='freq', yaxis='tsys', spw='17:4~123,19:4~123,21:4~123,23:4~123', iteraxis='antenna')
#---- WVR
wvrgcal(segsource=True, caltable=prefix + '.WVR', vis=prefix+'.ms', wvrspw=[4], toffset=0, wvrflag=[], statsource=PHCAL)
#---- Apply Tsys and WVR
tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable=prefix + 'T0.tsys', tsysChanTol=1)
for fields in list(set([BPCAL, PHCAL, TARGET])):
    applycal(vis=prefix+'.ms', field=fields, flagbackup=False, spw='25,27,29,31', interp='linear', gaintable=[prefix+'T0.tsys', prefix + '.WVR'], gainfield=fields, spwmap=[tsysmap, []], calwt=True)
os.system('rm -rf ' + prefix+'_Cal.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_Cal.ms', datacolumn='corrected', spw='25,27,29,31')
#-------- Phase Cal for bandpass
os.system('rm -rf P0*')
gaincal(vis=prefix+'_Cal.ms', caltable='P0', spw='*', field= BPCAL, scan='', selectdata=True, solint='int', refant=REFANT, calmode='p')
plotms('P0', xaxis='time', yaxis='phase', spw='*', iteraxis='antenna')
os.system('rm -rf B0*')
bandpass(vis = prefix + '_Cal.ms', caltable = 'B0', gaintable='P0', spw='*', field=  BPCAL, scan='', minblperant=5, minsnr=4, solint='inf', combine='scan', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True, spwmap=[0,1,2,3])
plotms('B0', xaxis='frequency', yaxis='amp', spw='*', iteraxis='antenna')
#-------- Phase Cal to align SPW 
os.system('rm -rf P1')
gaincal(vis=prefix+'_Cal.ms', caltable='P1', spw='*', gaintable = ['P0','B0'], field=BPCAL, selectdata=True, solint='inf,240ch', refant=REFANT, gaintype='G', calmode='p', minsnr=7, spwmap=[[1,1,1,1],[0,1,2,3]])
#-------- Phase Cal for all
os.system('rm -rf P2')
gaincal(vis=prefix+'_Cal.ms', caltable='P2', spw='*', gaintable = ['B0','P0','P1'], field='0,1', selectdata=True, solint='int', refant=REFANT, gaintype='G', combine='spw', calmode='p', minsnr=3, spwmap=[[0,1,2,3],[1,1,1,1],[0,1,2,3]])
plotms('P2', xaxis='time', yaxis='phase', spw='*', iteraxis='antenna', coloraxis='spw')
#-------- Flux cal
#Rscript ~/ALMA_Py3/polQuery.R -F314.0 -D2022-05-12T15:16:12 J0423-0120
#J0423-0120 2.91764981987018 -0.0562234647453836 -0.0620933289883079
IQUV = [2.91764981987018, -0.0562234647453836, -0.0620933289883079, 0.0]
setjy(vis=prefix + '_Cal.ms', field=BPCAL, spw='', standard='manual', fluxdensity=IQUV, spix = [-0.7,0], reffreq = '314.0GHz', usescratch=False)
os.system('rm -rf G0*')
gaincal(vis=prefix + '_Cal.ms', caltable = 'G0', spw ='', field=BPCAL + ',' + PHCAL, minsnr=3.0, solint='61s', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['P0', 'B0', 'P2'], calmode = 'a', gaintype='T', spwmap=[[0,1,2,3],[0,1,2,3],[0,0,0,0]], interp=['nearest', 'nearest','nearest'])
plotms('G0', xaxis = 'time', yaxis = 'amp', spw='*', iteraxis = 'antenna', coloraxis='corr')
fluxscale(vis=prefix + '_Cal.ms', caltable= 'G0', fluxtable='G0.flux', reference=BPCAL, transfer=PHCAL, refspwmap=[0,1,2,3])
#2022-05-26 14:45:29 INFO fluxscale	##### Begin Task: fluxscale          #####
#2022-05-26 14:45:29 INFO fluxscale	fluxscale( vis='uid___A002_Xf8d822_X9d6a_Cal.ms', caltable='G0', fluxtable='G0.flux', reference=['J0423-0120'], transfer=['J0241-0815'], listfile='', append=False, refspwmap=[0, 1, 2, 3], gainthreshold=-1.0, antenna='', timerange='', scan='', incremental=False, fitorder=1, display=False )
#2022-05-26 14:45:29 INFO fluxscale	****Using NEW VI2-driven calibrater tool****
#2022-05-26 14:45:29 INFO fluxscale	Opening MS: uid___A002_Xf8d822_X9d6a_Cal.ms for calibration.
#2022-05-26 14:45:29 INFO fluxscale	Initializing nominal selection to the whole MS.
#2022-05-26 14:45:29 INFO fluxscale	Beginning fluxscale--(MSSelection version)-------
#2022-05-26 14:45:29 INFO fluxscale	 Found reference field(s): J0423-0120
#2022-05-26 14:45:29 INFO fluxscale	 Found transfer field(s):  J0241-0815
#2022-05-26 14:45:30 INFO fluxscale	 Flux density for J0241-0815 in SpW=0 (freq=3.06409e+11 Hz) is: 0.401371 +/- 0.00596604 (SNR = 67.276, N = 40)
#2022-05-26 14:45:30 INFO fluxscale	 Flux density for J0241-0815 in SpW=1 (freq=3.08209e+11 Hz) is: 0.401411 +/- 0.00596048 (SNR = 67.3455, N = 40)
#2022-05-26 14:45:30 INFO fluxscale	 Flux density for J0241-0815 in SpW=2 (freq=3.18409e+11 Hz) is: 0.380952 +/- 0.00626688 (SNR = 60.7881, N = 40)
#2022-05-26 14:45:30 INFO fluxscale	 Flux density for J0241-0815 in SpW=3 (freq=3.20209e+11 Hz) is: 0.385951 +/- 0.00627873 (SNR = 61.4696, N = 40)
#2022-05-26 14:45:30 INFO fluxscale	 Fitted spectrum for J0241-0815 with fitorder=1: Flux density = 0.392372 +/- 0.00197175 (freq=313.25 GHz) spidx: a_1 (spectral index) =-1.14236 +/- 0.259284 covariance matrix for the fit:  covar(0,0)=6.10171e-05 covar(0,1)=0.000681988 covar(1,0)=0.000681988 covar(1,1)=0.861243
#2022-05-26 14:45:30 INFO fluxscale	Storing result in G0.flux
#2022-05-26 14:45:30 INFO fluxscale	Writing solutions to table: G0.flux
#2022-05-26 14:45:30 INFO fluxscale	Task fluxscale complete. Start time: 2022-05-26 10:45:29.276793 End time: 2022-05-26 10:45:29.841899
#2022-05-26 14:45:30 INFO fluxscale	##### End Task: fluxscale            #####
applycal(vis=prefix + '_Cal.ms', field='', flagbackup=False, spw='', interp=['nearest', 'nearest', 'nearest', 'nearest'], gaintable=['P0', 'B0', 'P2','G0.flux'], gainfield=[BPCAL, BPCAL, '',''], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]], calwt=True)
for sourceName in [BPCAL, PHCAL]:
    os.system('rm -rf ' + sourceName + '.ms')
    split(prefix + '_Cal.ms', spw='0,1,2,3', field=sourceName, outputvis=sourceName + '.ms', datacolumn='corrected')
#
#-------- CLEAN BPCAL
os.system('rm -rf ' + BPCAL + '.clean*') 
tclean(vis= BPCAL + '.ms', imagename=BPCAL + '.clean', field=BPCAL, cell=['0.1arcsec'], imsize=[128,128], spw='0,1,2,3', specmode='mfs', deconvolver='clark', weighting='briggs', robust=0.5, interactive=True, niter=1000, pbcor=False)
#-------- CLEAN NGC1052 continuum
os.system('rm -rf P3*') 
gaincal(vis=PHCAL+'.ms', caltable='P3', spw='*', selectdata=True, solint='int', refant=REFANT, gaintype='G', calmode='p', minsnr=3)
applycal(vis=PHCAL+'.ms', flagbackup=False, spw='*', interp=['nearest'], gaintable='P3', calwt=True)
plotms(PHCAL + '.ms', spw='0,1,2,3', xaxis='frequency', yaxis='amp', antenna='*&', avgtime='1e6', coloraxis='spw')
os.system('rm -rf ' + PHCAL + '.ContLSB.clean*') 
tclean(vis= PHCAL + '.ms', datacolumn='corrected', imagename=PHCAL + '.ContLSB.clean', field=PHCAL, cell=['0.1arcsec'], imsize=[256,256], spw='0,1', specmode='mfs', deconvolver='clark', weighting='natural', interactive=True, niter=1000, pbcor=False)
os.system('rm -rf ' + PHCAL + '.ContUSB.clean*') 
tclean(vis= PHCAL + '.ms', datacolumn='corrected', imagename=PHCAL + '.ContUSB.clean', field=PHCAL, cell=['0.1arcsec'], imsize=[256,256], spw='2:2~238,3:2~25;90~238', specmode='mfs', deconvolver='clark', weighting='natural', interactive=True, niter=1000, restoringbeam='', pbcor=False)
os.system('rm -rf ' + PHCAL + '.ContLine.SPW3.clean*') 
tclean(vis= PHCAL + '.ms', datacolumn='corrected', imagename=PHCAL + '.ContLine.SPW3.clean', field=PHCAL, cell=['0.1arcsec'], imsize=[128,128], spw='3', specmode='cube', start=0, nchan=240, width=1, outframe='LSRK', veltype='radio', restfreq='321.2256760GHz', deconvolver='clark', weighting='natural', interactive=True, niter=10000, pbcor=False)
os.system('rm -rf ' + PHCAL + '.ContLine.SPW2.clean*') 
tclean(vis= PHCAL + '.ms', datacolumn='corrected', imagename=PHCAL + '.ContLine.SPW2.clean', field=PHCAL, cell=['0.1arcsec'], imsize=[128,128], spw='2', specmode='cube', start=0, nchan=240, width=1, outframe='LSRK', veltype='radio', restfreq='321.2256760GHz', deconvolver='clark', weighting='natural', interactive=True, niter=10000, pbcor=False)
imview(PHCAL + '.ContLine.SPW2.clean.image')
#-------- UVCONTSUB
os.system('rm -rf ' + PHCAL + '*.contsub') 
uvcontsub(PHCAL + '.ms', fitspw='2:1~238,3:1~20;100~238', combine='spw', fitorder=1, spw='2,3')
os.system('rm -rf ' + PHCAL + '.Line.SPW3.clean*') 
tclean(vis= PHCAL + '.ms.contsub', datacolumn='corrected', imagename=PHCAL + '.Line.SPW3.clean', field=PHCAL, cell=['0.05arcsec'], imsize=[256,256], spw='1', specmode='cube', start='1225km/s', nchan=16, width='25km/s', outframe='LSRK', veltype='radio', restfreq='321.2256760GHz', deconvolver='hogbom', weighting='natural', interactive=True, restoringbeam='0.25arcsec', niter=1000, mosweight=True, pbcor=False, savemodel='modelcolumn')
#-------- Output maser and continuum visivilities
os.system('rm -rf ' + PHCAL + '_Maser.ms*') 
split(PHCAL + '.ms.contsub', outputvis=PHCAL + '_Maser.ms', spw='1:34~75', datacolumn  = 'data', width='42')
os.system('rm -rf ' + PHCAL + '_ContUSB.ms*') 
split(PHCAL + '.ms', outputvis=PHCAL + '_ContUSB.ms', spw='2:1~238', datacolumn  = 'corrected', width='238')
