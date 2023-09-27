from skimage.transform import iradon
import os
from os.path import expanduser, basename
import subprocess
import numpy as np
import sys
import h5py
import glob
import time

timestr = time.strftime("%Y%m%d%H%M%S")

def findNextPowerOf2(np2):
 np2=np2-1
 while np2&np2-1:
  np2=np2&np2-1
 return np2<<1

if len(sys.argv)!=2:
 print("Usage:\ anapy recon_peak_all_intg_only_v2_filter.py paramFN")
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
    if line.startswith('Ext'):
        ext = line.split()[1]
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
    #if line.startswith('OutFolder'):
        #OutFolder = line.split()[1]
    if line.startswith('ExtraPadForTomo'):
        extraPad = int(line.split()[1])
    if line.startswith('SeedFolder'):
        AnalyStem = line.split()[1]
    
if multipeak==1:
    rWidth=int((RMax-RMin)/2)
    
nFiles = endNr - startNr + 1
#fStem = RawFolder+FileStem
OutFolder=f'{AnalyStem}'+f'Output_{FileStem}_fileNr_{startNr}_{endNr}'
os.makedirs(OutFolder,exist_ok=True)
#UniqueId=f'RMin{RMin}_RMax{RMax}_RBinSize{rBinSize*100}%100_EtaMin{EtaMin}_EtaMax{EtaMax}_etaBinSize{etaBinSize*100}%100'

paramFN_split=paramFN.split('/')[-1]
updF = open(f'{OutFolder}/{paramFN_split}_intg_{timestr}', 'w')

paramContents = open(paramFN, 'r').readlines()
for line in paramContents:
    updF.write(line)
 
updF.write(f'timestr {timestr}\n')
updF.write(f'nFiles {nFiles}\n')
updF.write(f'OutFolder Output_{FileStem}_fileNr_{startNr}_{endNr}\n')
updF.write(f'OutFolder_Integration {OutFolder}'+f'/Integration\n')
updF.write(f'OutFolder_Peakfit {OutFolder}'+f'/Peakfit\n')
#updF.write(f'UniqueId {UniqueId}\n')

for rad in rads:
    updF.write(f'RadiusToFit {rad} {rWidth}\n')
for eta in etas:
    updF.write(f'EtaToFit {eta} {etaWidth}\n')
for Rcenter in Rcenters:
    updF.write(f'Rcenter {Rcenter}\n')
    
if multipeak==0:
   Rcenters=rads
updF.close()

#baseFN = basename(fStem)
#outfStem = f'{OutFolder}/{baseFN}'
#outfStem_integration=f'{OutFolder}'+f'/Integration/{baseFN}'
#outfStem_integration_Lineout_raw=f'{OutFolder}'+f'/Integration/Lineout_raw/{baseFN}'
#outfStem_integration_Lineout_wobg=f'{OutFolder}'+f'/Integration/Lineout_withoutbg/{baseFN}'
#outfStem_integration_intg=f'{OutFolder}'+f'/Integration/Integrated/{baseFN}'
#outfStem_integration_intg_1d=f'{OutFolder}'+f'/Integration/Integrated_1d/{baseFN}'
#outfStem_integration_REtaArea=f'{OutFolder}'+f'/Integration/REtaMap/{baseFN}'

#outfStem_peakfit=f'{OutFolder}'+f'/Peakfit/{baseFN}'


os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_mean',exist_ok=True)
os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_mean_withoutbg',exist_ok=True)
os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_median',exist_ok=True)
#os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_mean_2theta',exist_ok=True)
#os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_mean_withoutbg_2theta',exist_ok=True)
#os.makedirs(f'{OutFolder}'+f'/Integration/Lineout_median_2theta',exist_ok=True)



os.makedirs(f'{OutFolder}'+f'/Integration/Integrated',exist_ok=True)
os.makedirs(f'{OutFolder}'+f'/Integration/Integrated_1d',exist_ok=True)
os.makedirs(f'{OutFolder}'+f'/Integration/REtaMap',exist_ok=True)
#os.makedirs(f'{OutFolder}'+f'/Peakfit',exist_ok=True)

#############DO THE LINEOUTS#####

cmd1 = f'{expanduser("~/opt/MIDAS/DT/bin/DetectorMapper")} {OutFolder}/{paramFN_split}_intg_{timestr}'
subprocess.call(cmd1,shell=True)

cmd = f'{expanduser("~/opt/MIDAS/DT/binmultipeak/IntegratorOMP_lineout_v2_filter")} {OutFolder}/{paramFN_split}_intg_{timestr}'
subprocess.call(cmd,shell=True)
'''
for fnr in range(startNr,endNr+1):
 fn = f'{outfStem_integration}_{str(fnr).zfill(pad)}{ext}.LineOuts_withoutbg.bin'
 newname=f'{outfStem_integration_Lineout_wobg}_{str(fnr).zfill(pad)}{ext}.LineOuts_withoutbg.bin'
 os.rename(fn,newname)
 fn = f'{outfStem_integration}_{str(fnr).zfill(pad)}{ext}.LineOuts.bin'
 newname=f'{outfStem_integration_Lineout_raw}_{str(fnr).zfill(pad)}{ext}.LineOuts.bin'
 os.rename(fn,newname)
 fn = f'{outfStem_integration}_{str(fnr).zfill(pad)}{ext}.REtaAreaMap.csv'
 newname=f'{outfStem_integration_REtaArea}_{str(fnr).zfill(pad)}{ext}.REtaAreaMap.csv'
 os.rename(fn,newname)
 fn = f'{outfStem_integration}_{str(fnr).zfill(pad)}{ext}_integrated.bin'
 newname=f'{outfStem_integration_intg}_{str(fnr).zfill(pad)}{ext}_integrated.bin'
 os.rename(fn,newname)
 fn = f'{outfStem_integration}_{str(fnr).zfill(pad)}{ext}_integrated_1d.bin'
 newname=f'{outfStem_integration_intg_1d}_{str(fnr).zfill(pad)}{ext}_integrated_1d.bin'
 os.rename(fn,newname)
'''
