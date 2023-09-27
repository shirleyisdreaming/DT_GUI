
from skimage.transform import iradon
import os
from os.path import expanduser, basename
import subprocess
import numpy as np
import sys
import h5py
import glob
import time




def findNextPowerOf2(np2):
 np2=np2-1
 while np2&np2-1:
  np2=np2&np2-1
 return np2<<1

if len(sys.argv)!=2:
 print("Usage:\ anapy recon_peak_all_recon_only_v2_filter.py paramFN.upd")
paramFN=sys.argv[1]
rads=[]
etas=[]
Rcenters=[]
paramContents = open(paramFN, 'r').readlines()

for line in paramContents:
  
    if line.startswith('RawFolder'):
        RawFolder = line.split()[1]
    if line.startswith('FileStem'):
        FileStem = line.split()[1]
    #if line.startswith('darkFN'):
        #darkFN = line.split()[1]
    if line.startswith('Ext '):
        ext = str(line.split()[1])
    #if line.startswith('paramFN'):
        #paramFN = line.split()[1]
    if line.startswith('startNr'):
        startNr = int(line.split()[1])
    if line.startswith('endNr'):
        endNr = int(line.split()[1])
    if line.startswith('Padding'):
        pad = int(line.split()[1])
    if line.startswith('nFrames'):
        nFrames = int(line.split()[1])
    if line.startswith('startOme'):
        startOme = float(line.split()[1])
    if line.startswith('omeStep'):
        omeStep = float(line.split()[1])
    if line.startswith('BadRotation'):
        BadRotation = int(line.split()[1])
    if line.startswith('RMin'):
        RMin = float(line.split()[1])
    if line.startswith('RMax'):
        RMax = float(line.split()[1])
    #if line.startswith('radRange'):
        #radRange = line.split()[1]
    if line.startswith('rads'):
        rads_str = line.split()[1:]
        for rad_str in rads_str:
            if(rad_str.startswith("#")):
               break
            rads.append(int(rad_str))
            
    if line.startswith('Rwidth'):
        rWidth = line.split()[1]
    if line.startswith('RBinSize'):
        rBinSize = float(line.split()[1])
    #if line.startswith('etaRange'):
       # etaRange = line.split()[1]
    if line.startswith('EtaMin'):
        EtaMin = float(line.split()[1])
    if line.startswith('EtaMax'):
        EtaMax = float(line.split()[1])
    if line.startswith('etas'):
        etas_str = line.split()[1:]
        for eta_str in etas_str:
            etas.append(int(eta_str))
    if line.startswith('etaWidth'):
        etaWidth = line.split()[1]
    if line.startswith('EtaBinSize'):
        etaBinSize = float(line.split()[1])
    if line.startswith('filt'):
        filt = line.split()[1]
    if line.startswith('multipeak'):
        multipeak = int(line.split()[1])
    if line.startswith('Rcenters'):
        Rcenters_str = line.split()[1:]
        for Rcenter_str in Rcenters_str:
            Rcenters.append(int(Rcenter_str))
    #if line.startswith('fitWidth'):
       # fitWidth = line.split()[1]
    if line.startswith('numProcs'):
        numProcs = line.split()[1]
    
    if line.startswith('ExtraPadForTomo'):
        extraPad = int(line.split()[1])
    if line.startswith('SeedFolder'):
        AnalyStem = line.split()[1]
    if line.startswith('timestr'):
        timestr = line.split()[1]
    if line.startswith('Filter'):
        Filter = int(line.split()[1])
    if line.startswith('withoutremovebkgtomo'):
        worbkgtomo = int(line.split()[1])
    
if multipeak==1:
    rWidth=int((RMax-RMin)/2)
nFiles = endNr - startNr + 1

if multipeak==0:
   Rcenters=rads
   
OutFolder=f'{AnalyStem}'+f'Output_{FileStem}_fileNr_{startNr}_{endNr}'

outfStem = f'{OutFolder}/{FileStem}'
outfStem_integration=f'{OutFolder}'+f'/Integration/{FileStem}'
outfStem_integration_Lineout_raw=f'{OutFolder}'+f'/Integration/Lineout_mean/{FileStem}'
outfStem_integration_Lineout_wobg=f'{OutFolder}'+f'/Integration/Lineout_mean_withoutbg/{FileStem}'
outfStem_integration_Lineout_median=f'{OutFolder}'+f'/Integration/Lineout_median/{FileStem}'
outfStem_integration_intg=f'{OutFolder}'+f'/Integration/Integrated/{FileStem}'
outfStem_integration_intg_1d=f'{OutFolder}'+f'/Integration/Integrated_1d/{FileStem}'
outfStem_integration_REtaArea=f'{OutFolder}'+f'/Integration/REtaMap/{FileStem}'
outfStem_peakfit=f'{OutFolder}'+f'/Peakfit/{FileStem}'

if(worbkgtomo==1&Filter==1):
 OutString='raw_median'
if(worbkgtomo==1&Filter==0):
 OutString='raw_mean'
if(worbkgtomo==0&Filter==0):
 OutString='wobkg_mean'


outputs = ['RMEAN', 'MixFactor', 'SigmaG', 'SigmaL', 'MaxInt', 'BGFit',
           'BGSimple', 'MeanError', 'FitIntegratedIntensity', 'TotalIntensity', 'MaxIntensityObs','TotalIntensityBackgroundCorr']

##################Read Lineouts:
nEtas = len(etas)
nRads = len(rads)
noutput=len(outputs)
nRcenters=len(Rcenters)

fn = f'{outfStem_integration_Lineout_wobg}_{str(startNr).zfill(pad)}{ext}_LineOuts_mean_withoutbg_{timestr}.bin'
###########METHOD 1##############
#################################
##f = open(fn)
nEls=int(os.path.getsize(fn)/(8.0*nFrames))

nElsPerRad = int(nEls/nRads)

updF = open(paramFN, 'a')
updF.write(f'nElsTot {nEls}\n')
updF.write(f'nElsPerRad {nElsPerRad}\n')
updF.close()

if extraPad==1:
 reconSize=2*findNextPowerOf2(nFiles)
else:
 reconSize=findNextPowerOf2(nFiles)

sino_arr = np.zeros((nFiles,nFrames,nEtas,nElsPerRad*nRads))
for fnr in range(startNr,endNr+1):
 
 if(worbkgtomo==1&Filter==1):
  fn = f'{outfStem_integration_Lineout_median}_{str(fnr).zfill(pad)}{ext}.LineOuts_median_{timestr}.bin'
 if(worbkgtomo==1&Filter==0):
  fn = f'{outfStem_integration_Lineout_raw}_{str(fnr).zfill(pad)}{ext}.LineOuts_mean_{timestr}.bin'
 if(worbkgtomo==0&Filter==0):
  fn = f'{outfStem_integration_Lineout_wobg}_{str(fnr).zfill(pad)}{ext}_LineOuts_mean_withoutbg_{timestr}.bin'
 f = open(fn)
 #nEls = int(np.fromfile(f,dtype=np.ulonglong,count=(3))[0])
 if BadRotation == 1:
  if (fnr-startNr) % 2 == 1:
   data = np.fromfile(f,dtype=np.double,count=(nEls*nEtas*nFrames))
   data=np.reshape(data,(nFrames,nEtas,nEls))# Intensity values for the lineouts
   data=np.flip(data,0)
   data=np.reshape(data,(nFrames,nEtas,nElsPerRad*nRads))
  else:
   data = np.fromfile(f,dtype=np.double,count=(nEls*nEtas*nFrames)).reshape(nFrames,nEtas,nElsPerRad*nRads)
 else:
  data = np.fromfile(f,dtype=np.double,count=(nEls*nEtas*nFrames)).reshape(nFrames,nEtas,nElsPerRad*nRads)
 sino_arr[fnr-startNr] = data
sino_arr_flipped = np.transpose(sino_arr,(2,3,1,0)).astype('float32')
sino_arr_flipped.tofile(f'{OutFolder}/sinos_{OutString}_{timestr}.bin')

os.makedirs(f'{OutFolder}/Tomo_{OutString}',exist_ok=True)
fn = f'{OutFolder}/tomo_config_{OutString}_{timestr}.txt'
f = open(fn,'w')
f.write('dataFileName '+'/'+str(OutFolder)+f'/sinos_{OutString}_{timestr}.bin\n')
f.write('reconFileName '+'/'+str(OutFolder)+f'/Tomo_{OutString}/recon\n')
f.write('areSinos 1\n')
f.write('detXdim '+str(nFiles)+'\n')
f.write('detYdim '+str(nRads*nElsPerRad*nEtas)+'\n')
f.write('filter '+filt+'\n')
f.write('thetaRange '+str(startOme)+' '+str(startOme+(nFrames-1)*omeStep)+' '+str(omeStep)+'\n')
f.write('slicesToProcess -1\n')
f.write('shiftValues 0.000000 0.000000 0.500000\n')
f.write('ringRemovalCoefficient 1.0\n')
f.write('doLog 0\n')
f.write('ExtraPad '+str(extraPad)+'\n')
f.close()


cmdTomo = f'~/opt/MIDAS/TOMO/bin_debug/MIDAS_TOMO  {OutFolder}/tomo_config_{OutString}_{timestr}.txt {numProcs}'
subprocess.call(cmdTomo,shell=True)

recons = np.empty((nRads*nElsPerRad*nEtas,reconSize,reconSize))
for radN in range(nRads):
 for ElsPerRadN in range(nElsPerRad):
   for etaN in range(nEtas):
    fNr=(radN*nElsPerRad+ElsPerRadN)*nEtas+etaN
    os.makedirs(f'{OutFolder}/Tomo_{OutString}/Rad_{rads[radN]}/Eta_{etas[etaN]}',exist_ok=True)
    fn = f'/{OutFolder}/Tomo_{OutString}/recon_'+str(fNr).zfill(5)+'_001_p0000.0_'+str(reconSize)+'_'+str(reconSize)+'_float32.bin'
    f = open(fn)
    data = np.fromfile(f,dtype=np.float32,count=(reconSize*reconSize)).reshape((reconSize,reconSize))
    data[data<0]=0
    recons[fNr] = data
    newname=f'{OutFolder}/Tomo_{OutString}/Rad_{rads[radN]}/Eta_{etas[etaN]}/{FileStem}_tomo_Rad_{rads[radN]}_pm_{rWidth}_Eta_{etas[etaN]}_pm_{etaWidth}_RBinNr_{ElsPerRadN}_'+str(reconSize)+'x'+str(reconSize)+f'_float32_{timestr}.bin'
    
    data.astype(np.float32).tofile(newname)
    os.remove(fn)

recons_reshape = 1*np.transpose(recons.reshape((nRads*nElsPerRad,nEtas,reconSize,reconSize)))
recons_reshape = np.transpose(recons_reshape,(1,0,2,3))
recons_reshape.astype(np.double).tofile(f'/{OutFolder}/{OutString}_RawDataforPeakFit_{timestr}.bin')

updF = open(paramFN, 'a')
updF.write('RawDataPeakFN '+str(OutFolder)+f'/{OutString}_RawDataforPeakFit_{timestr}.bin\n')
updF.write('PeakFitResultFN '+str(OutFolder)+f'/{OutString}_PeakFitResult_{timestr}.bin\n')

updF.write('ReconSize '+str(reconSize)+'\n')
updF.close()



###peakfit give results of both methods###
print("Starting peak fitting")
subprocess.call('~/opt/MIDAS/DT/binmultipeak/PeakFit '+paramFN,shell=True)
print("Peak fitting finished")
if multipeak==0:
   fitResult = np.fromfile(f'{OutFolder}/{OutString}_PeakFitResult_{timestr}.bin',dtype=np.double,count=(reconSize*reconSize*nEtas*nRads*12)).reshape((reconSize,reconSize,nEtas,nRads,12)).transpose()
   for outputN in range(noutput):
     for radN in range(nRads):
         for etaN in range(nEtas):
             os.makedirs(f'{OutFolder}/Reconstruction/Rad_{rads[radN]}/Eta_{etas[etaN]}',exist_ok=True)
             fn = f'{OutFolder}/Reconstruction/Rad_{rads[radN]}/Eta_{etas[etaN]}/{FileStem}_FileNrs_{startNr}_{endNr}_{outputs[outputN]}_Rad_{rads[radN]}_pm_{rWidth}_Eta_{etas[etaN]}_pm_{etaWidth}_Rcenter_{Rcenters[radN]}_filter_{filt}_{reconSize}x{reconSize}_float64_{timestr}.bin'
             fitResult[outputN][radN][etaN].transpose().tofile(fn)
 
if multipeak==1:
    fitResult = np.fromfile(f'{OutFolder}/PeakFitResult_{timestr}.bin',dtype=np.double,count=(reconSize*reconSize*nEtas*nRads*nRcenters*12)).reshape((reconSize,reconSize,nEtas,nRads,nRcenters,12)).transpose()
    for outputN in range(noutput):
     for radN in range(nRads):
       for etaN in range(nEtas):
        for RcenterN in range(nRcenters):
            os.makedirs(f'{OutFolder}/Reconstruction/Rad_{rads[radN]}/Eta_{etas[etaN]}/Rcenter_{Rcenters[RcenterN]}',exist_ok=True)
            fn = f'{OutFolder}/Reconstruction/Rad_{rads[radN]}/Eta_{etas[etaN]}/Rcenter_{Rcenters[RcenterN]}/{FileStem}_FileNrs_{startNr}_{endNr}_{outputs[outputN]}_Rad_{rads[radN]}_pm_{rWidth}_Eta_{etas[etaN]}_pm_{etaWidth}_Rcenter_{Rcenters[RcenterN]}_filter_{filt}_{reconSize}x{reconSize}_float64_{timestr}.bin'
            fitResult[outputN][RcenterN][radN][etaN].transpose().tofile(fn)
 
'''

hf=h5py.File(f'{outfStem}'+'_DT_Recon.h5','w')

g1=hf.create_group('Integration')
g2=hf.create_group('Sino')
g3=hf.create_group('Tomo')
g4=hf.create_group('Reconstruction')

#data=[]
#for fnr in range(startNr,endNr+1):
    #data_now = np.fromfile(f'{outfStem}_{str(fnr).zfill(pad)}{ext}.REtaAreaMap.csv',dtype=np.float32)
    #data.append(data_now)
#g1.create_dataset('REtaAreaMap',data=data)

data=[]
data=np.fromfile(f'{OutFolder}/Integration/AreaFractionPixels.bin',dtype=np.float32)
g1.create_dataset('AreaFractionPixels',data=data)

nRBin=int((RMax-RMin)/rBinSize)
nEtaBin=int((EtaMax-EtaMin)/etaBinSize)
data=[]
for fnr in range(startNr,endNr+1):
    data_now = np.fromfile(f'{outfStem_integration_intg}_{str(fnr).zfill(pad)}{ext}_integrated.bin',dtype=np.double).reshape(nFrames,nRBin,nEtaBin)
    data.append(data_now)
g1.create_dataset('Integrated',data=data)

#data=[]
#for fnr in range(startNr,endNr+1):
    #data_now = np.fromfile(f'{outfStem_integration_intg_1d}_{str(fnr).zfill(pad)}{ext}_integrated_1d.bin',dtype=np.float32)
    #data.append(data_now)
#g1.create_dataset('Integrated_1d',data=data)

data=[]
for fnr in range(startNr,endNr+1):
    data_now = np.fromfile(f'{outfStem_integration_Lineout_raw}_{str(fnr).zfill(pad)}{ext}.LineOuts.bin',dtype=np.double).reshape(nFrames,nElsPerRad,nRads,nEtas)
    data.append(data_now)
g1.create_dataset('Lineout_raw',data=data)

data=[]
for fnr in range(startNr,endNr+1):
    data_now = np.fromfile(f'{outfStem_integration_Lineout_wobg}_{str(fnr).zfill(pad)}{ext}.LineOuts_withoutbg.bin',dtype=np.double).reshape(nFrames,nElsPerRad,nRads,nEtas)
    data.append(data_now)
g1.create_dataset('Lineout_BGcorrected',data=data)

data=[]
data=np.fromfile(f'{OutFolder}/sinos.bin',dtype=np.float32).reshape(nEtas,nRads,nElsPerRad,nFrames,nFiles)
g2.create_dataset('Sino',data=data)

data=[]
for rad in rads:
    for eta in etas:
        g3_subgroup=hf.create_group('Tomo/rad_'+f'{rad}'+'/'+'eta_'+f'{eta}')
        for ElsPerRadN in range(nElsPerRad):
            data_now=np.fromfile(f'{AnalyStem}/{OutFolder}/Tomo/Rad_{rad}/Eta_{eta}/recon_Rad_{rad}_pm_{rWidth}_Eta_{eta}_pm_{etaWidth}_RBinNr_{ElsPerRadN}_'+str(reconSize)+'x'+str(reconSize)+'_float32.bin',dtype=np.float32).reshape(reconSize,reconSize)
            data.append(data_now)
        g3_subgroup.create_dataset('rad_'+f'{rad}'+'_eta_'+f'{eta}',data=data)

data=[]
if multipeak==0:
   for rad in rads:
       for eta in etas:
           g4_subgroup=hf.create_group('Reconstruction/rad_'+f'{rad}'+'/'+'eta_'+f'{eta}')
           for output in outputs:
               data_now=np.fromfile(f'{AnalyStem}/{OutFolder}/Reconstruction/Rad_{rad}/Eta_{eta}/{baseFN}_FileNrs_{startNr}_{endNr}_{output}_Rad_{rad}_pm_{rWidth}_Eta_{eta}_pm_{etaWidth}_Rcenter_{rad}_filter_{filt}_{reconSize}x{reconSize}_float64.bin',dtype=np.double).reshape(reconSize,reconSize)
               data.append(data_now)
           g4_subgroup.create_dataset('rad_'+f'{rad}'+'_eta_'+f'{eta}',data=data)
  
data=[]
if multipeak==1:
   for rad in rads:
       for eta in etas:
           for Rcenter in Rcenters:
               g4_subgroup=hf.create_group('Reconstruction/rad_'+f'{rad}'+'/'+'eta_'+f'{eta}'+'/'+'Rcenter_'+f'{Rcenter}')
               for output in outputs:
                   data_now=np.fromfile(f'{AnalyStem}/{OutFolder}/Reconstruction/Rad_{rad}/Eta_{eta}/{baseFN}_FileNrs_{startNr}_{endNr}_{output}_Rad_{rad}_pm_{rWidth}_Eta_{etas}_pm_{etaWidth}_Rcenter_{Rcenter}_filter_{filt}_{reconSize}x{reconSize}_float64.bin',dtype=np.double).reshape(reconSize,reconSize)
                   data.append(data_now)
               g4_subgroup.create_dataset('rad_'+f'{rad}'+'_eta_'+f'{eta}'+'_Rcenter_'+f'{Rcenter}',data=data)

'''
