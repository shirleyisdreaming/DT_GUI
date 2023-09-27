import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import sys
import tkinter as Tk
from tkinter import ttk
import numpy as np
import time
import matplotlib.pyplot as plt
import os
import tkinter.filedialog as tkFileDialog
import math
import scipy
import scipy.ndimage as ndimage
import subprocess
from subprocess import call
from matplotlib import gridspec
from os.path import expanduser, basename
from textwrap import fill

def _quit():
    root.quit()
    root.destroy()
 
def getfn_integration():
 file_integration=int(file_integration_var.get())
 return fstem_intg+'/'+OutFolder_intg+'/Integration/Integrated/'+samplename_intg+'_'+str(file_integration).zfill(6)+ext_intg+'_integrated_'+timestr_intg+'.bin'
 
def getfn_integration_1d():
 file_integration=int(file_integration_var.get())
 return fstem_intg+'/'+OutFolder_intg+'/Integration/Integrated_1d/'+samplename_intg+'_'+str(file_integration).zfill(6)+ext_intg+'_integrated_1d_'+timestr_intg+'.bin'
  
def getfn_lineout_raw():
 file_lineout=int(file_lineout_var.get())
 return fstem+'/'+OutFolder+'/Integration/Lineout_mean/'+samplename+'_'+str(file_lineout).zfill(6)+ext+'_LineOuts_mean_'+timestr+'.bin'
 
def getfn_lineout_withoutbg():
 file_lineout=int(file_lineout_var.get())
 return fstem+'/'+OutFolder+'/Integration/Lineout_mean_withoutbg/'+samplename+'_'+str(file_lineout).zfill(6)+ext+'_LineOuts_mean_withoutbg_'+timestr+'.bin'
 
def getfn_sino():
 return fstem+'/'+OutFolder+'/sinos_'+timestr+'.bin'
 
def getfn_tomo():
 global elementNr_tomo,rad_tomo,eta_tomo
 elementNr_tomo = int(elementNrvar_tomo.get())
 return fstem+'/'+OutFolder+'/Tomo/Rad_'+f'{rad_tomo}'+'/Eta_'+f'{eta_tomo}/'+samplename+'_tomo_'+'Rad_'+str(rad_tomo)+'_pm_'+str(Rwidth)+'_Eta_'+str(eta_tomo)+'_pm_'+str(etaWidth)+'_RBinNr_'+str(elementNr_tomo)+'_'+str(reconSize)+'x'+str(reconSize)+'_float32_'+timestr+'.bin'

def getfn_recon():
 global Resulttype_recon,rad_recon,eta_recon,Rcenter_recon
 if multipeak==0:
    Rcenter_recon=rad_recon
    return fstem+'/'+OutFolder+'/Reconstruction/Rad_'+f'{rad_recon}'+'/Eta_'+f'{eta_recon}'+'/'+samplename+'_FileNrs_'+str(StartNr)+'_'+str(EndNr)+'_'+str(Resulttype_recon)+'_Rad_'+str(rad_recon)+'_pm_'+str(Rwidth)+'_Eta_'+str(eta_recon)+'_pm_'+str(etaWidth)+'_Rcenter_'+str(Rcenter_recon)+'_filter_'+str(filterNr)+'_'+str(reconSize)+'x'+str(reconSize)+'_float64_'+timestr+'.bin'
 if multipeak==1:
    return fstem+'/'+OutFolder+'/Reconstruction/Rad_'+f'{rad_recon}'+'/Eta_'+f'{eta_recon}'+'/Rcenter_'+f'{Rcenter_recon}'+'/'+samplename+'_FileNrs_'+str(StartNr)+'_'+str(EndNr)+'_'+str(Resulttype_recon)+'_Rad_'+str(rad_recon)+'_pm_'+str(Rwidth)+'_Eta_'+str(eta_recon)+'_pm_'+str(etaWidth)+'_Rcenter_'+str(Rcenter_recon)+'_filter_'+str(filterNr)+'_'+str(reconSize)+'x'+str(reconSize)+'_float64_'+timestr+'.bin'


def getData_integration():
    global frame_integration
    frame_integration=int(frame_integration_var.get())
    fn=getfn_integration()
    print("Reading Integrated File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.double,count=(nRBin_intg*nEtaBin_intg*nFrames_intg))
    f.close()
    data = np.reshape(data,(nFrames_intg,nRBin_intg,nEtaBin_intg))
    data=data[frame_integration,:,:]
    data=np.transpose(data,(1,0))
    #data=np.reshape(data,(nRBin_intg,nEtaBin_intg))
    return data
    
def getData_integration_1d():
    global frame_integration
    frame_integration=int(frame_integration_var.get())
    fn=getfn_integration_1d()
    print("Reading Integrated_1d File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,sep=" ",dtype=np.float64)
    f.close()
    data = np.reshape(data,(nFrames_intg,nRBin_intg,3))
    data=data[frame_integration,:,2]
    return data
    
def getData_lineout_raw():
    global eta_lineout,rad_lineout,frame_lineout
    frame_lineout=int(frame_lineout_var.get())
    fn=getfn_lineout_raw()
    print("Reading Raw Lineout File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.double,count=(nElsTot*nEtaFit*nFrames))
    f.close()
    data = np.reshape(data,(nFrames,nEtaFit,nElsTot))
    eta_lineout_Nr=0
    rad_lineout_Nr=0
    for i in range(len(Eta)):
      if Eta[i]==eta_lineout:
       eta_lineout_Nr=i
    for j in range(len(Rad)):
      if Rad[j]==rad_lineout:
       rad_lineout_Nr=j
    data=data[frame_lineout,eta_lineout_Nr,nElsPerRad*(rad_lineout_Nr):nElsPerRad*(rad_lineout_Nr+1)]
    return data
    
def getData_lineout_withoutbg():
    global eta_lineout,rad_lineout,frame_lineout
    frame_lineout=int(frame_lineout_var.get())
    fn=getfn_lineout_withoutbg()
    print("Reading Lineout without Background File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.double,count=(nElsTot*nEtaFit*nFrames))
    f.close()
    data = np.reshape(data,(nFrames,nElsTot,nEtaFit))
    eta_lineout_Nr=0
    rad_lineout_Nr=0
    for i in range(len(Eta)):
      if Eta[i]==eta_lineout:
       eta_lineout_Nr=i
    for j in range(len(Rad)):
      if Rad[j]==rad_lineout:
       rad_lineout_Nr=j
    data=data[frame_lineout,nElsPerRad*(rad_lineout_Nr):nElsPerRad*(rad_lineout_Nr+1),eta_lineout_Nr]
    return data
    
def getData_sino():
    global eta_sino,rad_sino,elementNr_sino
    elementNr_sino=int(elementNrvar_sino.get())
    fn=getfn_sino()
    print("Reading Sino File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.float32,count=(nEtaFit*nElsTot*nFiles*nFrames))
    f.close()
    data = np.reshape(data,(nEtaFit,nElsPerRad*nRads,nFrames,nFiles))
    eta_sino_Nr=0
    rad_sino_Nr=0
    for i in range(len(Eta)):
      if Eta[i]==eta_sino:
       eta_sino_Nr=i
    for j in range(len(Rad)):
      if Rad[j]==rad_sino:
       rad_sino_Nr=j
    data=data[eta_sino_Nr,(nElsPerRad*(rad_sino_Nr)):(nElsPerRad*(rad_sino_Nr+1)),:,:]
    data=data[(elementNr_sino):(elementNr_sino+1),:,:]
    return data
    
def getData_tomo():
    fn=getfn_tomo()
    print("Reading Tomo File: " + fn+" "+str(elementNr_tomo))
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.float32,count=(reconSize*reconSize))
    f.close()
    data = np.reshape(data,(reconSize,reconSize))
    return data

def getData_tomo_all():
    data_all=[]
    for elementNr_now in range(elementNr):
     fn=fstem+'/'+OutFolder+'/Tomo/Rad_'+f'{rad_tomo}'+'/Eta_'+f'{eta_tomo}/'+samplename+'_tomo_'+'Rad_'+str(rad_tomo)+'_pm_'+str(Rwidth)+'_Eta_'+str(eta_tomo)+'_pm_'+str(etaWidth)+'_RBinNr_'+str(elementNr_now)+'_'+str(reconSize)+'x'+str(reconSize)+'_float32_'+timestr+'.bin'
     f = open(fn,'rb')
     data = np.fromfile(f,dtype=np.float32,count=(reconSize*reconSize))
     f.close()
     #data = np.reshape(data,(reconSize,reconSize))
     data_all.append(data)
    data_all=np.reshape(data_all,(elementNr,reconSize,reconSize))
    return data_all
    
def getData_recon():
    fn=getfn_recon()
    print("Reading Reconstruction File: " + fn)
    f = open(fn,'rb')
    data = np.fromfile(f,dtype=np.double,count=(reconSize*reconSize))
    f.close()
    data = np.reshape(data,(reconSize,reconSize))
    return data

def Loadplot_integration():
    global data_integration,data_integration_1d
    threshold = float(thresholdvar_integration.get())
    upperthreshold = float(maxthresholdvar_integration.get())
    file_integration = int(file_integration_var.get())
    frame_integration = int(frame_integration_var.get())
    data_integration = getData_integration()
    data_integration_1d=getData_integration_1d()

    integration.clear()
    integration_1d.clear()
    integration_table.clear()
    
    pl_intg=integration.imshow(data_integration,cmap=plt.get_cmap('hot'),interpolation='nearest',clim=(threshold,upperthreshold),extent=[Rmin_intg,Rmax_intg,Etamin_intg,Etamax_intg])
    
   #pl_intg=integration.imshow(data_integration,cmap=plt.get_cmap('hot'),interpolation='nearest',clim=(threshold,upperthreshold),extent=[math.degrees(math.atan(Rmin_intg*172.0/1071096.0)),math.degrees(math.atan(Rmin_intg*172.0/1071096.0)),Etamin_intg,Etamax_intg])
   
    integration.title.set_text("Integration")
    integration.set_xlabel("2Theta")
    integration.set_ylabel("Eta Value")
    integration.set_aspect('auto')

    cax = figur_integration.add_axes([0.9,0.75,0.01,0.2])
    c=figur_integration.colorbar(pl_intg,cax=cax)
    c.set_label('Intensity Counts')
   
    x=np.linspace(math.degrees(math.atan(Rmin_intg*172.0/1071096.0)),math.degrees(math.atan(Rmax_intg*172.0/1071096.0)),nRBin_intg)
    integration_1d.axis(xmin=math.degrees(math.atan(Rmin_intg*172.0/1071096.0)),xmax=math.degrees(math.atan(Rmax_intg*172.0/1071096.0)))
    integration_1d.scatter(x,data_integration_1d)
    integration_1d.title.set_text("Integration_1d")
    integration_1d.set_xlabel("2Theta")
    integration_1d.set_ylabel("Intensity Counts")
    
    row_labels=['paramFN_Intg','sample','fileNr','Start Omega','Omega Step','Omega Now','RMin~RMax','RBinSize','nRBin','EtaMin~EtaMax','EtaBinSize','nEtaBin']
    #textwrap.wrap(filePath_param_intg_read,width=6)
    temp_w=fill(filePath_param_intg_read,width=70)
    table_vals=[[temp_w],[samplename_intg],[file_integration],[startOme_intg],[omeStep_intg],[startOme_intg+frame_integration*omeStep_intg],[f'{Rmin_intg}'+'~'+f'{Rmax_intg}'],[RBinsize_intg],[nRBin_intg],[f'{Etamin_intg}'+'~'+f'{Etamax_intg}'],[EtaBinsize_intg],[nEtaBin_intg]]
    
    tb=integration_table.table(cellText=table_vals,rowLabels=row_labels,loc='center',cellLoc='center')
    cellDict=tb.get_celld()
    cellDict[(0,0)].set_height(.2)
    
    integration_table.set_title("Parameters for Integration")
    tb.auto_set_column_width([0,1])
    #tb.scale(0.5,1.0)
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    integration_table.axis('off')
    
    canvas_integration.draw()

def incr_file_plotupdater_integration():
    file_integration = int(file_integration_var.get())
    file_integration += 1
    file_integration_var.set(str(file_integration))
    Loadplot_integration()

def decr_file_plotupdater_integration():
    file_integration = int(file_integration_var.get())
    file_integration -= 1
    file_integration_var.set(str(file_integration))
    Loadplot_integration()
    
def incr_frame_plotupdater_integration():
    frame_integration = int(frame_integration_var.get())
    frame_integration += 1
    frame_integration_var.set(str(frame_integration))
    Loadplot_integration()
    
def decr_frame_plotupdater_integration():
    frame_integration = int(frame_integration_var.get())
    frame_integration -= 1
    frame_integration_var.set(str(frame_integration))
    Loadplot_integration()

def Loadplot_lineout():
    global data_lineout_raw,data_lineout_withoutbg
    file_lineout = int(file_lineout_var.get())
    data_lineout_raw = getData_lineout_raw()
    data_lineout_withoutbg = getData_lineout_withoutbg()
    lineout.clear()
    
    x=np.arange(rad_lineout-Rwidth-RBinSize,rad_lineout+Rwidth+RBinSize,RBinSize)
    lineout.scatter(x,data_lineout_raw)
    lineout.scatter(x,data_lineout_withoutbg)
    
    lineout.legend(["raw data","Background Removal"])
    lineout.set_xlabel("radii")
    lineout.set_ylabel("Intensity")
    lineout.title.set_text("Lineout Display")
    
    temp_w=fill(filePath_param,width=70)
    if multipeak==0:
     row_labels=['paramFN_lineout','sample','fileNr','Start Omega','Omega Step','Omega Now','Radii','Radii Range','RBinSize','nRBin','Eta','Eta Range','EtaBinSize']
    
     table_vals=[[temp_w],[samplename],[file_lineout],[startOme],[omeStep],[startOme+frame_lineout*omeStep],[rad_lineout],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[eta_lineout],['±'+f'{etaWidth}'],[EtaBinSize]]
    if multipeak==1:
     row_labels=['paramFN_lineout','sample','fileNr','Start Omega','Omega Step','Omega Now','multipeak','Radii','Radii Range','Peak Centers','RBinSize','nRBin','Eta','Eta Range','EtaBinSize']
     table_vals=[[temp_w],[samplename],[file_lineout],[startOme],[omeStep],[startOme+frame_lineout*omeStep],['Yes'],[rad_lineout],['±'+f'{Rwidth}'],[Rcenters],[RBinSize],[nElsPerRad],[eta_lineout],['±'+f'{etaWidth}'],[EtaBinSize]]
    
    tb=lineout_table.table(cellText=table_vals,rowLabels=row_labels,loc='center',cellLoc='center')
    cellDict=tb.get_celld()
    cellDict[(0,0)].set_height(.2)
    
    lineout_table.set_title("Parameters for lineout")
    tb.auto_set_column_width([0,1])
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    lineout_table.axis('off')

    canvas_lineout.draw()

def incr_plotupdater_file_lineout():
    file_lineout = int(file_lineout_var.get())
    file_lineout += 1
    file_lineout_var.set(str(file_lineout))
    Loadplot_lineout()

def decr_plotupdater_file_lineout():
    file_lineout = int(file_lineout_var.get())
    file_lineout -= 1
    file_lineout_var.set(str(file_lineout))
    Loadplot_lineout()

def incr_plotupdater_frame_lineout():
    frame_lineout = int(frame_lineout_var.get())
    frame_lineout += 1
    frame_lineout_var.set(str(frame_lineout))
    Loadplot_lineout()

def decr_plotupdater_frame_lineout():
    frame_lineout = int(frame_lineout_var.get())
    frame_lineout -= 1
    frame_lineout_var.set(str(frame_lineout))
    Loadplot_lineout()

def Loadplot_sino():
    global data_sino
    threshold = float(thresholdvar_sino.get())
    upperthreshold = float(maxthresholdvar_sino.get())
    elementNr_sino=int(elementNrvar_sino.get())
    data_sino = getData_sino()
    sino.imshow(data_sino.reshape(nFrames,nFiles).transpose(1,0),cmap=plt.get_cmap('bone'),interpolation='nearest',clim=(threshold,upperthreshold),extent=[startOme,startOme+nFrames*omeStep,StartNr,StartNr+nFiles-1])
    #sino.imshow(data_sino[(elementNr_sino):(elementNr_sino+1),:,:].reshape(nFrames,nFiles).transpose(1,0),cmap=plt.get_cmap('bone'),interpolation='nearest',clim=(threshold,upperthreshold),extent=[startOme,startOme+nFrames*omeStep,StartNr,StartNr+nFiles-1])
    sino.set_xlabel("Omega value")
    sino.set_ylabel("File Nr")
    sino.title.set_text("Sinogram Display")
    temp_w=fill(filePath_param,width=70)
    if multipeak==0:
     row_labels=['paramFN_sino','sample','file Range','Omega Range','Radii','Radii Range','RBinSize','nRBin','Eta','Eta Range','EtaBinSize','Radii now']
     table_vals=[[temp_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],[rad_sino],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[eta_sino],['±'+f'{etaWidth}'],[EtaBinSize],[rad_sino-Rwidth+RBinSize*elementNr_sino]]
    if multipeak==1:
     row_labels=['paramFN_sino','sample','file Range','Omega Range','Fit Multipeaks','Radii','Radii Range','RBinSize','nRBin','Rcenters','Eta','Eta Range','EtaBinSize','Radii now']
     table_vals=[[temp_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],['Yes'],[rad_sino],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[Rcenters],[eta_sino],['±'+f'{etaWidth}'],[EtaBinSize],[rad_sino-Rwidth+RBinSize*elementNr_sino]]
     
    tb=sino_table.table(cellText=table_vals,rowLabels=row_labels,loc='center',cellLoc='center')
    
    cellDict=tb.get_celld()
    cellDict[(0,0)].set_height(.2)
    sino_table.set_title("Parameters for sinogram")
    #tb.auto_set_column_width(1,2)
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    sino_table.axis('off')
    
    canvas_sino.draw()

def incr_plotupdater_sino():
    elementNr_sino = int(elementNrvar_sino.get())
    elementNr_sino += 1
    elementNrvar_sino.set(str(elementNr_sino))
    Loadplot_sino()

def decr_plotupdater_sino():
    elementNr_sino = int(elementNrvar_sino.get())
    elementNr_sino -= 1
    elementNrvar_sino.set(str(elementNr_sino))
    Loadplot_sino()
 
 
def onclick(event):
    if event.dblclick:
     tomo_line.clear()
     xdata=int(event.xdata)
     ydata=int(event.ydata)
     data_tomo_line=data_tomo_all[:,xdata,ydata]
     #x=np.arange(0,elementNr,1)
     x=np.arange(rad_tomo-Rwidth-RBinSize,rad_tomo+Rwidth+RBinSize,RBinSize)
     tomo_line.scatter(x,data_tomo_line)
     tomo_line.set_title("Lineout Contributes to Voxel: x=%d y=%d"%(xdata,ydata))
     tomo_line.set_xlabel("Radii")
     tomo_line.set_ylabel("Intensity")
     canvas_tomo.draw()
     
def Loadplot_tomo():
    global data_tomo,elementNr_tomo,data_tomo_all
    threshold = float(thresholdvar_tomo.get())
    upperthreshold = float(maxthresholdvar_tomo.get())
    elementNr_tomo = int(elementNrvar_tomo.get())
    data_tomo = getData_tomo()
    data_tomo_all=getData_tomo_all()
    refreshPlot = 1
    tomo.imshow(data_tomo,cmap=plt.get_cmap('bone'),interpolation='nearest',clim=(threshold,upperthreshold))
    tomo.title.set_text("Tomography Display")
    tomo_line.text(0.1,0.1,'Please double click the Tomo Figure \nto probe the voxel \nfor checking the lineout')
    canvas_tomo.mpl_connect('button_press_event', onclick)
    temp_w=fill(filePath_param,width=70)
    if multipeak==0:
     row_labels=['paramFN_tomo','sample','file Range','Omega Range','Rad','Rad Range','RBinSize','nRBin','Eta','Eta Range','EtaBinSize','reconstruction size','Radii now']
    
     table_vals=[[temp_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],[rad_tomo],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[eta_tomo],['±'+f'{etaWidth}'],[EtaBinSize],[reconSize],[rad_tomo-Rwidth+RBinSize*elementNr_tomo]]
    if multipeak==1:
     row_labels=['paramFN_tomo','sample','file Range','Omega Range','Fit Multipeaks','Rad','Rad Range','RBinSize','nRBin','Rcenters','Eta','Eta Range','EtaBinSize','reconstruction size','Radii now']
     table_vals=[[temp_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],['Yes'],[rad_tomo],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[Rcenters],[eta_tomo],['±'+f'{etaWidth}'],[EtaBinSize],[reconSize],[rad_tomo-Rwidth+RBinSize*elementNr_tomo]]
     
    tb=tomo_table.table(cellText=table_vals,rowLabels=row_labels,loc='center',cellLoc='center')
    cellDict=tb.get_celld()
    cellDict[(0,0)].set_height(.2)
    tomo_table.set_title("Parameters for Tomography")
    tb.auto_set_column_width([0,1])
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    tomo_table.axis('off')

    canvas_tomo.draw()
   
def incr_plotupdater_tomo():
    elementNr_tomo = int(elementNrvar_tomo.get())
    elementNr_tomo += 1
    elementNrvar_tomo.set(str(elementNr_tomo))
    Loadplot_tomo()
    
def decr_plotupdater_tomo():
    elementNr_tomo = int(elementNrvar_tomo.get())
    elementNr_tomo -= 1
    elementNrvar_tomo.set(str(elementNr_tomo))
    Loadplot_tomo()

def Loadplot_recon():
    global refreshPlot
    global data_recon
    threshold = float(thresholdvar_m1_recon.get())
    upperthreshold = float(maxthresholdvar_m1_recon.get())
    data_recon = getData_recon()
    refreshPlot = 1
    m1_recon.imshow(data_recon,cmap=plt.get_cmap('bone'),interpolation='nearest',clim=(threshold,upperthreshold))
    m1_recon.title.set_text("Reconstruction Display")
    temp_w=fill(filePath_param,width=70)
    if multipeak==0:
     row_labels=['paramFN_recon','sample','file Range','Omega Range','Radii','Radii Range','RBinSize','nRBin','Eta','Eta Range','EtaBinSize','Rcenter choosed now']
    
     table_vals=[[tem_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],[rad_recon],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[eta_recon],['±'+f'{etaWidth}'],[EtaBinSize],[Rcenter_recon]]
    if multipeak==1:
     row_labels=['paramFN_recon','sample','file Range','Omega Range','Fit Multipeaks','Radii','Radii Range','RBinSize','nRBin','Rcenters','Eta','Eta Range','EtaBinSize','Rcenter choosed now']
     table_vals=[[temp_w],[samplename],[f'{StartNr}'+'~'+f'{StartNr+nFiles-1}'],[f'{startOme}'+'~'+f'{startOme+omeStep*nFrames}'],['Yes'],[rad_recon],['±'+f'{Rwidth}'],[RBinSize],[nElsPerRad],[Rcenters],[eta_recon],['±'+f'{etaWidth}'],[EtaBinSize],[Rcenter_recon]]
     
    tb=recon_table.table(cellText=table_vals,rowLabels=row_labels,loc='center',cellLoc='center')
    cellDict=tb.get_celld()
    cellDict[(0,0)].set_height(.2)
    recon_table.set_title("Parameters for Reconstruction")
    #tb.auto_set_column_width(1,2)
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    recon_table.axis('off')

    
    canvas_recon.draw()
    
def readParams():
    global filePath_param
    global StartNr,EndNr,samplename,nFrames,Rad,Rwidth,Eta,etaWidth,nFiles,filterNr,elementNr,fileStem,nElsTot,fstem,nEtaFit,nRads,nElsPerRad,reconSize,OutFolder,nRBin,nEtaBin,multipeak,RBinSize,EtaBinSize,Rcenters,startOme,omeStep,timestr,Lsd,ext
    print("Load parameters from "+filePath_param+" to this DT workflow")
    f= open(filePath_param,'r')
    paramContents=f.readlines()
    Eta=[]
    Rad=[]
    FileNr_value=[]
    Rcenters=[]
    nRads=0
    nEtaFit=0
    nRcenter=0
    for line in paramContents:
        if line == '\n':
            continue
        if line[0] == '#':
            continue
                    
        if 'startNr' == line.split()[0]:
            StartNr = int(line.split()[1])
        if 'endNr' == line.split()[0]:
            EndNr = int(line.split()[1])
        if 'FileStem' == line.split()[0]:
            samplename = str((line.split()[1]))
        if 'timestr'==line.split()[0]:
            timestr=str(line.split()[1])
        if 'nFrames' == line.split()[0]:
            nFrames = int(line.split()[1])
        if 'RadiusToFit'==line.split()[0]:
            Rad.append(int(line.split()[1]))
            Rwidth=int(line.split()[-1])
            nRads=nRads+1
        if 'Rcenters'==line.split()[0]:
            rcs_str = line.split()[1:]
            for rc_str in rcs_str:
                Rcenters.append(int(rc_str))
                nRcenter=nRcenter+1
        if 'EtaToFit'==line.split()[0]:
            Eta.append(int(line.split()[1]))
            etaWidth=int(line.split()[-1])
            nEtaFit=nEtaFit+1
        if 'filt'==line.split()[0]:
            filterNr=int(line.split()[1])
        if 'nElsPerRad'==line.split()[0]:
            elementNr=int(line.split()[1])
            nElsPerRad=int(line.split()[1])
        if 'ReconSize'==line.split()[0]:
            reconSize=int(line.split()[1])
        if 'nElsTot'==line.split()[0]:
            nElsTot=int(line.split()[1])
        if 'SeedFolder'==line.split()[0]:
            fstem=str(line.split()[1])
        if 'OutFolder'==line.split()[0]:
            OutFolder=str(line.split()[1])
        if 'nRBin'==line.split()[0]:
            nRBin=int(line.split()[1])
        if 'nEtaBin'==line.split()[0]:
            nEtaBin=int(line.split()[1])
        if 'multipeak'==line.split()[0]:
            multipeak=int(line.split()[1])
        if 'RBinSize'==line.split()[0]:
            RBinSize=float(line.split()[1])
        if 'EtaBinSize'==line.split()[0]:
            EtaBinSize=float(line.split()[1])
        if 'startOme'==line.split()[0]:
            startOme=float(line.split()[1])
        if 'omeStep'==line.split()[0]:
            omeStep=float(line.split()[1])
        if 'Lsd'==line.split()[0]:
            Lsd=float(line.split()[1])
        if 'Ext' == line.split()[0]:
            ext = str(line.split()[1])
            
            
    nFiles=EndNr-StartNr+1
    FileNr_value=np.linspace(StartNr,EndNr,1)
    cbox_rad_sino.config(value=Rad)
    cbox_eta_sino.config(value=Eta)
    cbox_rad_lineout.config(value=Rad)
    cbox_eta_lineout.config(value=Eta)
 
def readParams_intg():
    global filePath_param_intg_read
    global StartNr_intg,EndNr_intg,samplename_intg,nFrames_intg,nFiles_intg,fileStem_intg,nElsTot_intg,fstem_intg,OutFolder_intg,RBinsize_intg,EtaBinsize_intg,Rmin_intg,Rmax_intg,Etamin_intg,Etamax_intg,nRBin_intg,nEtaBin_intg,startOme_intg,omeStep_intg,timestr_intg,Lsd_intg,ext_intg
    print("Load parameters from "+filePath_param_intg_read+" for checking integration result")
    f= open(filePath_param_intg_read,'r')
    paramContents=f.readlines()
    for line in paramContents:
        if line == '\n':
            continue
        if line[0] == '#':
            continue
        if 'startNr' == line.split()[0]:
            StartNr_intg = int(line.split()[1])
        if 'endNr' == line.split()[0]:
            EndNr_intg = int(line.split()[1])
        if 'FileStem' == line.split()[0]:
            samplename_intg = str(line.split()[1])
        if 'nFrames' == line.split()[0]:
            nFrames_intg = int(line.split()[1])
        if 'SeedFolder'==line.split()[0]:
            fstem_intg=str(line.split()[1])
        if 'OutFolder'==line.split()[0]:
            OutFolder_intg=str(line.split()[1])
        if 'timestr'==line.split()[0]:
            timestr_intg=str(line.split()[1])
        if 'RMin'==line.split()[0]:
            Rmin_intg=int(line.split()[1])
        if 'RMax'==line.split()[0]:
            Rmax_intg=int(line.split()[1])
        if 'EtaMin'==line.split()[0]:
            Etamin_intg=int(line.split()[1])
        if 'EtaMax'==line.split()[0]:
            Etamax_intg=int(line.split()[1])
        if 'RBinSize'==line.split()[0]:
            RBinsize_intg=float(line.split()[1])
        if 'EtaBinSize'==line.split()[0]:
            EtaBinsize_intg=float(line.split()[1])
        if 'startOme'==line.split()[0]:
            startOme_intg=float(line.split()[1])
        if 'omeStep'==line.split()[0]:
            omeStep_intg=float(line.split()[1])
        if 'Lsd'==line.split()[0]:
            Lsd_intg=float(line.split()[1])
        if 'Ext' == line.split()[0]:
            ext_intg = str(line.split()[1])
            
    nRBin_intg=int((Rmax_intg-Rmin_intg)/RBinsize_intg)
    nEtaBin_intg=int((Etamax_intg-Etamin_intg)/EtaBinsize_intg)
    
    
def folderSelector_parameter_intg():
   global filePath_param_intg
   filePath_param_intg = Tk.StringVar()
   filePath_param_intg=Tk.filedialog.askopenfilename(title=u'choose file')
   print("Choosed "+filePath_param_intg+" for Integration")
   filePath_param_intg_en.grid(row=26,column=2,sticky=Tk.W)
   filePath_param_intg_en.delete(0,"end")
   filePath_param_intg_en.insert(0,filePath_param_intg)

def folderSelector_parameter_intg_read():
   global filePath_param_intg_read
   filePath_param_intg_read = Tk.StringVar()
   filePath_param_intg_read=Tk.filedialog.askopenfilename(title=u'choose file')
   print("Choosed "+filePath_param_intg_read+" for checking Integration")
   filePath_param_intg_read_en.grid(row=40,column=2,sticky=Tk.W)
   filePath_param_intg_read_en.delete(0,"end")
   filePath_param_intg_read_en.insert(0,filePath_param_intg_read)
   
def folderSelector_parameter_recon():
   global filePath_param_recon
   filePath_param_recon = Tk.StringVar()
   filePath_param_recon=Tk.filedialog.askopenfilename(title=u'choose file')
   print("Choosed "+filePath_param_recon+" for Reconstruction")
   filePath_param_recon_en.grid(row=50,column=2,sticky=Tk.W)
   filePath_param_recon_en.delete(0,"end")
   filePath_param_recon_en.insert(0,filePath_param_recon)
   
   
def folderSelector_parameter_wf():
   global filePath_param
   filePath_param = Tk.StringVar()
   filePath_param=Tk.filedialog.askopenfilename(title=u'choose file')
   print("Choosed "+filePath_param+" for Checking Reconstruction Results")
   filePath_param_en.grid(row=56,column=2,sticky=Tk.W)
   filePath_param_en.delete(0,"end")
   filePath_param_en.insert(0,filePath_param)

def func_lineout_rad(event):
    global rad_lineout
    rad_lineout=float(rad_lineout_var.get())
    
def func_lineout_eta(event):
    global eta_lineout
    eta_lineout=float(eta_lineout_var.get())

def func_sino_rad(event):
    global rad_sino
    rad_sino=float(rad_sino_var.get())
    
def func_sino_eta(event):
    global eta_sino
    eta_sino=float(eta_sino_var.get())
   
def func_tomo_rad(event):
    global rad_tomo
    rad_tomo=int(rad_tomo_var.get())

def func_tomo_eta(event):
    global eta_tomo
    eta_tomo=int(eta_tomo_var.get())

def func_recon_rstype(event):
    global Resulttype_recon
    Resulttype_recon=str(resultstyle_recon_var.get())

def func_recon_rad(event):
    global rad_recon
    rad_recon=int(rad_recon_var.get())

def func_recon_eta(event):
    global eta_recon
    eta_recon=int(eta_recon_var.get())
    
def func_recon_rcenter(event):
    global Rcenter_recon
    Rcenter_recon=int(Rcenter_recon_var.get())
    
def auto_threshold_sino():
    threshold_sino=np.min(data_sino)
    upperthreshold=np.max(data_sino)
    thresholdvar_sino.set(str(threshold_sino))
    maxthresholdvar_sino.set(str(upperthreshold))
    
def auto_threshold_tomo():
    threshold_tomo=np.min(data_tomo)
    upperthreshold=np.max(data_tomo)
    thresholdvar_tomo.set(str(threshold_tomo))
    maxthresholdvar_tomo.set(str(upperthreshold))

def auto_threshold_recon():
    threshold_recon=np.min(data_recon)
    upperthreshold=np.max(data_recon)
    thresholdvar_m1_recon.set(str(threshold_recon))
    maxthresholdvar_m1_recon.set(str(upperthreshold))


def show_integration_progress():
    global p_intg,label_intg
    cmd_intg='/home/beams12/S1IDUSER/opt/midasconda/bin/python'+' '+'/home/beams12/S1IDUSER/opt/MIDAS/DT/recon_peak_all_intg_only_v2_filter.py'+' '+str(filePath_param_intg)
    p_intg=subprocess.Popen("exec "+cmd_intg,shell=True)
    label_intg.grid_remove()
    progress_integration='Running: '+f'{filePath_param_intg}'
    label_intg=Tk.Label(master=firstRowFrame, text=progress_integration, font=("Helvetica",12))
    label_intg.grid(row=30,column=2,padx=10)
    root.update()
    a=p_intg.wait()
    if a>1:
       label_intg.grid_remove()
       progress_integration='Error. Please check the parameter file: '+f'{filePath_param_intg}'
       label_intg=Tk.Label(master=firstRowFrame, text=progress_integration, font=("Helvetica",12))
       label_intg.grid(row=30,column=2,padx=10)
       root.update()
    if a==0:
       label_intg.grid_remove()
       progress_integration='Finished: '+f'{filePath_param_intg}'
       label_intg=Tk.Label(master=firstRowFrame, text=progress_integration, font=("Helvetica",12))
       label_intg.grid(row=30,column=2,padx=10)
       root.update()
    
def stop_integration():
    p_intg.kill()
    print("Integration Progress stopped")
    label_intg.grid_remove()
    progress_integration='Stopped: '+f'{filePath_param_intg}'
    label_intg=Tk.Label(master=firstRowFrame, text=progress_integration, font=("Helvetica",12))
    label_intg.grid(row=30,column=2,padx=10)
    root.update()

def show_recon_progress():
    global p,label_recon
    cmd_recon='/home/beams12/S1IDUSER/opt/midasconda/bin/python'+' '+'/home/beams12/S1IDUSER/opt/MIDAS/DT/recon_peak_all_recon_only_v2_filter_2.py'+' '+str(filePath_param_recon)
    p=subprocess.Popen("exec "+cmd_recon,shell=True)
    label_recon.grid_remove()
    progress_Reconstruction='Running '+f'{filePath_param_recon}'
    label_recon=Tk.Label(master=firstRowFrame, text=progress_Reconstruction, font=("Helvetica",12))
    label_recon.grid(row=51,column=2,padx=10)
    root.update()
    a=p.wait()
    if a>1:
       label_recon.grid_remove()
       progress_Reconstruction='Error. Please check the parameter file: '+f'{filePath_param_recon}'
       label_recon=Tk.Label(master=firstRowFrame, text=progress_Reconstruction, font=("Helvetica",12))
       label_recon.grid(row=51,column=2,padx=10)
       print('Error. Please check the parameter file.')
       root.update()
    if a==0:
       label_recon.grid_remove()
       progress_Reconstruction='Finished: '+f'{filePath_param_recon}'
       label_recon=Tk.Label(master=firstRowFrame, text=progress_Reconstruction, font=("Helvetica",12))
       label_recon.grid(row=51,column=2,padx=10)
       root.update()


def stop_progress():
    p.kill()
    print("Reconstruction Progress stopped")
    label_recon.grid_remove()
    progress_Reconstruction='Stopped: '+f'{filePath_param_recon}'
    label_recon=Tk.Label(master=firstRowFrame, text=progress_Reconstruction, font=("Helvetica",12))
    label_recon.grid(row=51,column=2,padx=10)
    root.update()

def createNewWindow_integration():
     global integration,canvas_integration,integration_1d,integration_table,figur_integration,integration_cbar
     integrationWindow= Tk.Toplevel(root)
     integrationWindow.wm_title(" Check Integration: "+ f'{samplename_intg}' )
     figur_integration = Figure(figsize=(8,8),dpi=100)
     canvas_integration = FigureCanvasTkAgg(figur_integration,master=integrationWindow)
     
     spec=gridspec.GridSpec(ncols=1, nrows=3)
     
     integration_1d = figur_integration.add_subplot(spec[1])
     integration = figur_integration.add_subplot(spec[0])
     #integration_1d = figur_integration.add_subplot(spec[1],sharex=integration)
     integration_table=figur_integration.add_subplot(spec[2])
     
     
     figur_integration.subplots_adjust(top=0.95,bottom=0.05,hspace=0.4)
 
     figrowspan_integration = 10
     figcolspan_integration = 10

     refreshPlot = 0

     canvas_integration.get_tk_widget().grid(row=0,column=0,columnspan=figcolspan_integration,rowspan=figrowspan_integration,sticky=Tk.W+Tk.E+Tk.N+Tk.S)
     toolbar_frame_integration = Tk.Frame(integrationWindow)
     toolbar_frame_integration.grid(row=figrowspan_integration+4,column=0,columnspan=10,sticky=Tk.W)
     toolbar_integration = NavigationToolbar2Tk( canvas_integration, toolbar_frame_integration )
     toolbar_integration.update()

     threshold_integration = 0
     thresholdvar_integration.set(str(threshold_integration))
     maxthresholdvar_integration.set(str(100))

     firstRowFrame_integration = Tk.Frame(integrationWindow)
     firstRowFrame_integration.grid(row=figrowspan_integration+6,column=1,rowspan=10,columnspan=10,sticky=Tk.W)
  
     file_integration=StartNr_intg
     file_integration_var.set(str(file_integration))
     frame_integration=0
     frame_integration_var.set(str(frame_integration))
     Tk.Label(master=firstRowFrame_integration,text="FileNr",font=("Helvetica",12)).grid(row=1,column=2,padx=(10,0))
     Tk.Label(master=firstRowFrame_integration,text="FileNr:(%d~%d)"%(StartNr_intg,EndNr_intg),font=("Helvetica",12)).grid(row=2,column=2,padx=(10,0))
     Tk.Entry(master=firstRowFrame_integration,textvariable=file_integration_var,width=6).grid(row=1,column=3)
     Tk.Button(master=firstRowFrame_integration,text='+',command=incr_file_plotupdater_integration,font=("Helvetica",12)).grid(row=1,column=4)
     Tk.Button(master=firstRowFrame_integration,text='-',command=decr_file_plotupdater_integration,font=("Helvetica",12)).grid(row=1,column=5)
     
     Tk.Label(master=firstRowFrame_integration,text="frameNr",font=("Helvetica",12)).grid(row=1,column=6,padx=(10,0))
     Tk.Label(master=firstRowFrame_integration,text="frameNr:(0~%d)"%(nFrames_intg),font=("Helvetica",12)).grid(row=2,column=6,padx=(10,0))
     Tk.Entry(master=firstRowFrame_integration,textvariable=frame_integration_var,width=6).grid(row=1,column=7)
     Tk.Button(master=firstRowFrame_integration,text='+',command=incr_frame_plotupdater_integration,font=("Helvetica",12)).grid(row=1,column=8)
     Tk.Button(master=firstRowFrame_integration,text='-',command=decr_frame_plotupdater_integration,font=("Helvetica",12)).grid(row=1,column=9)
     
     Tk.Button(master=firstRowFrame_integration,text="Load Integration & Integration_1d",command=Loadplot_integration,font=("Helvetica",14)).grid(row=1,column=10,padx=(20,10))

     Tk.Label(master=firstRowFrame_integration,text='MinThresh',font=("Helvetica",12)).grid(row=1,column=11,padx=(10,0))
     Tk.Entry(master=firstRowFrame_integration,textvariable=thresholdvar_integration,width=5).grid(row=1,column=12,padx=(0,10))
     Tk.Label(master=firstRowFrame_integration,text='MaxThresh',font=("Helvetica",12)).grid(row=1,column=13,padx=(10,0))
     Tk.Entry(master=firstRowFrame_integration,textvariable=maxthresholdvar_integration,width=5).grid(row=1,column=14,padx=(0,10))
     Tk.Button(master=firstRowFrame_integration,text='UpdThresh',font=("Helvetica",12),command=Loadplot_integration).grid(row=1,column=15,padx=(10,0))


def createNewWindow_lineout():
     global lineout,canvas_lineout,lineout_table
     lineoutWindow= Tk.Toplevel(root)
     lineoutWindow.wm_title(" Check Lineout "+f'{samplename}')
     figur_lineout = Figure(figsize=(10,10),dpi=100)
     canvas_lineout = FigureCanvasTkAgg(figur_lineout,master=lineoutWindow)
     spec=gridspec.GridSpec(ncols=1, nrows=2)

     lineout = figur_lineout.add_subplot(spec[0])
     lineout_table=figur_lineout.add_subplot(spec[1])
     figur_lineout.subplots_adjust(left=0.2,hspace=0.5)
     
     figrowspan_lineout = 10
     figcolspan_lineout = 10

     refreshPlot = 0

     canvas_lineout.get_tk_widget().grid(row=0,column=0,columnspan=figcolspan_lineout,rowspan=figrowspan_lineout,sticky=Tk.W+Tk.E+Tk.N+Tk.S)
     toolbar_frame_lineout = Tk.Frame(lineoutWindow)
     toolbar_frame_lineout.grid(row=figrowspan_lineout+4,column=0,columnspan=10,sticky=Tk.W)
     toolbar_lineout = NavigationToolbar2Tk( canvas_lineout, toolbar_frame_lineout )
     toolbar_lineout.update()

     firstRowFrame_lineout = Tk.Frame(lineoutWindow)
     firstRowFrame_lineout.grid(row=figrowspan_lineout+6,column=1,rowspan=10,columnspan=10,sticky=Tk.W)
  
     file_lineout=StartNr
     file_lineout_var.set(str(file_lineout))
     Tk.Label(master=firstRowFrame_lineout,text="FileNr",font=("Helvetica",12)).grid(row=1,column=2,padx=(10,0))
     Tk.Label(master=firstRowFrame_lineout,text="FileNr:(%d~%d)"%(StartNr,EndNr),font=("Helvetica",12)).grid(row=2,column=2,padx=(10,0))
     
     Tk.Entry(master=firstRowFrame_lineout,textvariable=file_lineout_var,width=6).grid(row=1,column=3)
     Tk.Button(master=firstRowFrame_lineout,text='+',command=incr_plotupdater_file_lineout,font=("Helvetica",12)).grid(row=1,column=4)
     Tk.Button(master=firstRowFrame_lineout,text='-',command=decr_plotupdater_file_lineout,font=("Helvetica",12)).grid(row=1,column=5)

     frame_lineout=0
     frame_lineout_var.set(str(frame_lineout))
     Tk.Label(master=firstRowFrame_lineout,text="frameNr",font=("Helvetica",12)).grid(row=1,column=6,padx=(10,0))
     Tk.Label(master=firstRowFrame_lineout,text="frameNr:(0~%d)"%(nFrames),font=("Helvetica",12)).grid(row=2,column=6,padx=(10,0))
     Tk.Entry(master=firstRowFrame_lineout,textvariable=frame_lineout_var,width=6).grid(row=1,column=7)
     Tk.Button(master=firstRowFrame_lineout,text='+',command=incr_plotupdater_frame_lineout,font=("Helvetica",12)).grid(row=1,column=8)
     Tk.Button(master=firstRowFrame_lineout,text='-',command=decr_plotupdater_frame_lineout,font=("Helvetica",12)).grid(row=1,column=9)

     Tk.Label(master=firstRowFrame_lineout,text="Rad Value",font=("Helvetica",12)).grid(row=1,column=10,padx=(20,0))
     cbox_rad_lineout = ttk.Combobox(master=firstRowFrame_lineout,textvariable=rad_lineout_var,value=Rad)
     cbox_rad_lineout.grid(row=1,column=11)
     cbox_rad_lineout.bind("<<ComboboxSelected>>",func_lineout_rad)
     
     Tk.Label(master=firstRowFrame_lineout,text="Eta Value",font=("Helvetica",12)).grid(row=1,column=12,padx=(20,0))
     cbox_eta_lineout = ttk.Combobox(master=firstRowFrame_lineout,textvariable=eta_lineout_var,value=Eta)
     cbox_eta_lineout.grid(row=1,column=13)
     cbox_eta_lineout.bind("<<ComboboxSelected>>",func_lineout_eta)

     Tk.Button(master=firstRowFrame_lineout,text="Load Lineout",command=Loadplot_lineout,font=("Helvetica",14)).grid(row=1,column=14,padx=(20,10))

 
def createNewWindow_sino():
    global sino,canvas_sino,sino_table
    sinoWindow= Tk.Toplevel(root)
    sinoWindow.wm_title(" Check Sinogram "+f'{samplename}')
    figur_sino = Figure(figsize=(9,7),dpi=100)
    canvas_sino = FigureCanvasTkAgg(figur_sino,master=sinoWindow)
    
    spec=gridspec.GridSpec(ncols=1, nrows=2)

    sino = figur_sino.add_subplot(spec[0])
    sino_table=figur_sino.add_subplot(spec[1])
    figur_sino.subplots_adjust(left=0.2,hspace=0.5)
    
    figrowspan_sino = 10
    figcolspan_sino = 10
    refreshPlot = 0
    threshold_sino = 0
    thresholdvar_sino.set(str(threshold_sino))
    maxthresholdvar_sino.set(str(1))
    
    
    toolbar_frame_sino = Tk.Frame(sinoWindow)
    toolbar_frame_sino.grid(row=figrowspan_sino+4,column=0,columnspan=10,sticky=Tk.W)
    toolbar_sino = NavigationToolbar2Tk( canvas_sino, toolbar_frame_sino )
    toolbar_sino.update()

    firstRowFrame_sino = Tk.Frame(sinoWindow)
    firstRowFrame_sino.grid(row=figrowspan_sino+6,column=1)
    
    Tk.Label(master=firstRowFrame_sino,text="Rad Value",font=("Helvetica",12)).grid(row=1,column=2,padx=(20,0))
    cbox_rad_sino = ttk.Combobox(master=firstRowFrame_sino,textvariable=rad_sino_var,value=Rad)
    cbox_rad_sino.grid(row=1,column=3)
    cbox_rad_sino.bind("<<ComboboxSelected>>",func_sino_rad)
    
    Tk.Label(master=firstRowFrame_sino,text="Eta Value",font=("Helvetica",12)).grid(row=1,column=4,padx=(10,0))
    cbox_eta_sino = ttk.Combobox(master=firstRowFrame_sino,textvariable=eta_sino_var,value=Eta)
    cbox_eta_sino.grid(row=1,column=5)
    cbox_eta_sino.bind("<<ComboboxSelected>>",func_sino_eta)
  
    elementNrvar_sino.set(str(elementNr_sino))
    Tk.Label(master=firstRowFrame_sino,text="RBinNr",font=("Helvetica",12)).grid(row=1,column=6,padx=(20,0))
    Tk.Label(master=firstRowFrame_sino,text="RBinNr:(0~%d)"%(nElsPerRad),font=("Helvetica",12)).grid(row=2,column=6,padx=(20,0))
    Tk.Entry(master=firstRowFrame_sino,textvariable=elementNrvar_sino,width=6).grid(row=1,column=7)
    Tk.Button(master=firstRowFrame_sino,text='+',command=incr_plotupdater_sino,font=("Helvetica",12)).grid(row=1,column=8,padx=(0,0))
    Tk.Button(master=firstRowFrame_sino,text='-',command=decr_plotupdater_sino,font=("Helvetica",12)).grid(row=1,column=9,padx=(0,0))

    Tk.Label(master=firstRowFrame_sino,text='MinThresh',font=("Helvetica",12)).grid(row=1,column=10,padx=(10,0))
    Tk.Entry(master=firstRowFrame_sino,textvariable=thresholdvar_sino,width=5).grid(row=1,column=11,padx=(0,10))
    Tk.Label(master=firstRowFrame_sino,text='MaxThresh',font=("Helvetica",12)).grid(row=1,column=12,padx=(10,0))
    Tk.Entry(master=firstRowFrame_sino,textvariable=maxthresholdvar_sino,width=5).grid(row=1,column=13,padx=(0,10))
    Tk.Button(master=firstRowFrame_sino,text='UpdThresh',font=("Helvetica",12),command=Loadplot_sino).grid(row=1,column=14,padx=(10,0))
    
    Tk.Button(master=firstRowFrame_sino,text='Auto Threshold',command=auto_threshold_sino,font=("Helvetica",12)).grid(row=1,column=15,padx=(0,0))
    Tk.Button(master=firstRowFrame_sino,text="Load Sino",command=Loadplot_sino,font=("Helvetica",14)).grid(row=1,column=16,padx=(10,0))
    
    canvas_sino.get_tk_widget().grid(row=0,column=0,columnspan=figcolspan_sino,rowspan=figrowspan_sino,sticky=Tk.W+Tk.E+Tk.N+Tk.S)

           
def createNewWindow_tomo():
     global tomo,canvas_tomo,elementNr_tomo,tomo_table,tomo_line
     tomoWindow= Tk.Toplevel(root)
     tomoWindow.wm_title(" Check Tomography "+f'{samplename}')
     figur_tomo = Figure(figsize=(8,8),dpi=100)
     canvas_tomo = FigureCanvasTkAgg(figur_tomo,master=tomoWindow)
     spec=gridspec.GridSpec(ncols=2, nrows=2)
     
     tomo = figur_tomo.add_subplot(spec[0])
     tomo_line=figur_tomo.add_subplot(spec[1])
     tomo_table = figur_tomo.add_subplot(spec[2])
     figur_tomo.subplots_adjust(left=0.3,hspace=0.5)
     
     figrowspan_tomo = 10
     figcolspan_tomo = 10
     elementNr_tomo=0
     threshold_tomo = 0
     thresholdvar_tomo.set(str(threshold_tomo))
     maxthresholdvar_tomo.set(str(1))
     
     refreshPlot = 0
     
     toolbar_frame_tomo = Tk.Frame(tomoWindow)
     toolbar_frame_tomo.grid(row=figrowspan_tomo+4,column=0,columnspan=10,sticky=Tk.W)
     toolbar_tomo = NavigationToolbar2Tk( canvas_tomo, toolbar_frame_tomo )
     toolbar_tomo.update()

     firstRowFrame_tomo = Tk.Frame(tomoWindow)
     firstRowFrame_tomo.grid(row=figrowspan_tomo+6,column=1)

     Tk.Label(master=firstRowFrame_tomo,text="Eta Value",font=("Helvetica",12)).grid(row=1,column=2,padx=(10,0))
     cbox_Eta_tomo = ttk.Combobox(master=firstRowFrame_tomo,textvariable=eta_tomo_var,value=Eta)
     cbox_Eta_tomo.grid(row=1,column=3,padx=(0,10))
     cbox_Eta_tomo.bind("<<ComboboxSelected>>",func_tomo_eta)
     
     Tk.Label(master=firstRowFrame_tomo,text="Rad Value",font=("Helvetica",12)).grid(row=1,column=4,padx=(10,0))
     cbox_Rad_tomo = ttk.Combobox(master=firstRowFrame_tomo,textvariable=rad_tomo_var,value=Rad)
     cbox_Rad_tomo.grid(row=1,column=5,padx=(0,10))
     cbox_Rad_tomo.bind("<<ComboboxSelected>>",func_tomo_rad)
     
     elementNrvar_tomo.set(str(elementNr_tomo))
     Tk.Label(master=firstRowFrame_tomo,text="RBinNr",font=("Helvetica",12)).grid(row=1,column=6,padx=(10,0))
     Tk.Label(master=firstRowFrame_tomo,text="RBinNr:(0~%d)"%(nElsPerRad),font=("Helvetica",12)).grid(row=2,column=6,padx=(20,0))
     Tk.Entry(master=firstRowFrame_tomo,textvariable=elementNrvar_tomo,width=6).grid(row=1,column=7)
     Tk.Button(master=firstRowFrame_tomo,text='+',command=incr_plotupdater_tomo,font=("Helvetica",12)).grid(row=1,column=8)
     Tk.Button(master=firstRowFrame_tomo,text='-',command=decr_plotupdater_tomo,font=("Helvetica",12)).grid(row=1,column=9)
     
     Tk.Label(master=firstRowFrame_tomo,text='MinThresh',font=("Helvetica",12)).grid(row=1,column=10,padx=(20,0))
     Tk.Entry(master=firstRowFrame_tomo,textvariable=thresholdvar_tomo,width=6).grid(row=1,column=11)
     Tk.Label(master=firstRowFrame_tomo,text='MaxThresh',font=("Helvetica",12)).grid(row=1,column=12,padx=(20,0))
     Tk.Entry(master=firstRowFrame_tomo,textvariable=maxthresholdvar_tomo,width=6).grid(row=1,column=13)
     Tk.Button(master=firstRowFrame_tomo,text='UpdThresh',font=("Helvetica",12),command=Loadplot_tomo).grid(row=1,column=14,padx=(10,0))

     Tk.Button(master=firstRowFrame_tomo,text='Auto Threshold',command=auto_threshold_tomo,font=("Helvetica",12)).grid(row=1,column=15,padx=(0,0))
     Tk.Button(master=firstRowFrame_tomo,text="Load Tomo",command=Loadplot_tomo,font=("Helvetica",14)).grid(row=1,column=16,padx=(20,0))
     canvas_tomo.get_tk_widget().grid(row=0,column=0,columnspan=figcolspan_tomo,rowspan=figrowspan_tomo,sticky=Tk.W+Tk.E+Tk.N+Tk.S)
    
def createNewWindow_recon():
     global m1_recon,canvas_recon,recon_table
     reconWindow= Tk.Toplevel(root)
     reconWindow.wm_title(" Reconstruction Result "+f'{samplename}')
     figur_recon = Figure(figsize=(8,8),dpi=100)
     canvas_recon = FigureCanvasTkAgg(figur_recon,master=reconWindow)
     spec=gridspec.GridSpec(ncols=1, nrows=2)
     
     m1_recon = figur_recon.add_subplot(spec[0])
     recon_table = figur_recon.add_subplot(spec[1])
     figur_recon.subplots_adjust(left=0.2,hspace=0.5)
     
     figrowspan_recon = 10
     figcolspan_recon = 10

     threshold_recon = 0
     thresholdvar_m1_recon.set(str(threshold_recon))
     maxthresholdvar_m1_recon.set(str(1))

     Resultvalue=["RMEAN", "MixFactor", "SigmaG", "SigmaL", "MaxInt", "BGFit",
           "BGSimple", "MeanError", "FitIntegratedIntensity", "TotalIntensity", "TotalIntensityBackgroundCorr", "MaxIntensityObs"]

     refreshPlot = 0

     canvas_recon.get_tk_widget().grid(row=0,column=0,columnspan=figcolspan_recon,rowspan=figrowspan_recon,sticky=Tk.W+Tk.E+Tk.N+Tk.S)
     toolbar_frame_recon = Tk.Frame(reconWindow)
     toolbar_frame_recon.grid(row=figrowspan_recon+4,column=0,columnspan=10,sticky=Tk.W)
     toolbar_recon = NavigationToolbar2Tk( canvas_recon, toolbar_frame_recon )
     toolbar_recon.update()

     firstRowFrame_recon = Tk.Frame(reconWindow)
     firstRowFrame_recon.grid(row=figrowspan_recon+6,column=1,sticky=Tk.W)

     
     Tk.Label(master=firstRowFrame_recon,text="Eta Value",font=("Helvetica",12)).grid(row=1,column=2,padx=(10,0))
     cbox_Eta_recon = ttk.Combobox(master=firstRowFrame_recon,textvariable=eta_recon_var,value=Eta)
     cbox_Eta_recon.grid(row=1,column=3,padx=(0,10))
     cbox_Eta_recon.bind("<<ComboboxSelected>>",func_recon_eta)
     
     Tk.Label(master=firstRowFrame_recon,text="Rad Value",font=("Helvetica",12)).grid(row=1,column=4,padx=(10,0))
     cbox_Rad_recon = ttk.Combobox(master=firstRowFrame_recon,textvariable=rad_recon_var,value=Rad)
     cbox_Rad_recon.grid(row=1,column=5,padx=(0,10))
     cbox_Rad_recon.bind("<<ComboboxSelected>>",func_recon_rad)
     
     if(multipeak==1):
      Tk.Label(master=firstRowFrame_recon,text="Rcenter Value",font=("Helvetica",12)).grid(row=1,column=6,padx=(10,0))
      cbox_Rcenter_recon = ttk.Combobox(master=firstRowFrame_recon,textvariable=Rcenter_recon_var,value=Rcenters)
      cbox_Rcenter_recon.grid(row=1,column=7,padx=(0,10))
      cbox_Rcenter_recon.bind("<<ComboboxSelected>>",func_recon_rcenter)
      
     Tk.Label(master=firstRowFrame_recon,text="Fitting Type",font=("Helvetica",12)).grid(row=1,column=8,padx=(10,0))
     cbox_recon = ttk.Combobox(master=firstRowFrame_recon,textvariable=resultstyle_recon_var,value=Resultvalue)
     cbox_recon.grid(row=1,column=9,padx=(0,10))
     cbox_recon.bind("<<ComboboxSelected>>",func_recon_rstype)
     
     Tk.Label(master=firstRowFrame_recon,text='MinThresh',font=("Helvetica",12)).grid(row=1,column=10,padx=(10,0))
     Tk.Entry(master=firstRowFrame_recon,textvariable=thresholdvar_m1_recon,width=6).grid(row=1,column=11,padx=(0,10))
     Tk.Label(master=firstRowFrame_recon,text='MaxThresh',font=("Helvetica",12)).grid(row=1,column=12,padx=(10,0))
     Tk.Entry(master=firstRowFrame_recon,textvariable=maxthresholdvar_m1_recon,width=6).grid(row=1,column=13,padx=(0,10))
     Tk.Button(master=firstRowFrame_recon,text='UpdThresh',font=("Helvetica",12),command=Loadplot_recon).grid(row=1,column=14,padx=(10,10))
     
     Tk.Button(master=firstRowFrame_recon,text='Auto Threshold',command=auto_threshold_recon,font=("Helvetica",12)).grid(row=1,column=15,padx=(0,0))
     
     Tk.Button(master=firstRowFrame_recon,text='Load Recon',font=("Helvetica",14),command=Loadplot_recon).grid(row=1,column=16)

def linkToff_aysm():
    cmd1='/home/beams12/S1IDUSER/opt/midasconda/bin/python'+' '+'/home/beams12/S1IDUSER/opt/MIDAS/gui/ff_asym_link_v2.py'+' '+str(filePath_param_intg)
    subprocess.call(cmd1,shell=True)

def open_parameterfile_intg():
    global filePath_param_intg
    cmd='geany'+' '+str(filePath_param_intg)
    subprocess.call(cmd,shell=True)

def open_parameterfile_recon():
    global filePath_param_recon
    cmd='geany'+' '+str(filePath_param_recon)
    subprocess.call(cmd,shell=True)

    
root = Tk.Tk()
root.wm_title(" DT Reconstruction Workflow")
figrowspan = 1
figcolspan = 1

paramFN = 'filepara.txt'
Rad=[]
Eta=[]
FileNr_value=[]

refreshPlot = 0
progress_integration='No Task'
progress_Reconstruction='No Task'

thresholdvar_integration = Tk.StringVar()
maxthresholdvar_integration= Tk.StringVar()

thresholdvar_sino = Tk.StringVar()
maxthresholdvar_sino= Tk.StringVar()
elementNr_sino=0
elementNrvar_sino = Tk.StringVar()
eta_sino_var=Tk.StringVar()
rad_sino_var=Tk.StringVar()
cbox_rad_sino=ttk.Combobox()
cbox_eta_sino=ttk.Combobox()

thresholdvar_tomo = Tk.StringVar()
maxthresholdvar_tomo= Tk.StringVar()
elementNr_tomo=0
elementNrvar_tomo = Tk.StringVar()

maxthresholdvar_m1_recon= Tk.StringVar()
thresholdvar_m1_recon = Tk.StringVar()

file_integration_var=Tk.StringVar()
frame_integration_var=Tk.StringVar()

file_lineout=0
file_lineout_var=Tk.StringVar()
frame_lineout=0
frame_lineout_var=Tk.StringVar()
rad_lineout_var=Tk.StringVar()
eta_lineout_var=Tk.StringVar()
cbox_rad_lineout=ttk.Combobox()
cbox_eta_lineout=ttk.Combobox()

eta_recon_var=Tk.StringVar()
rad_recon_var=Tk.StringVar()
eta_tomo_var=Tk.StringVar()
rad_tomo_var=Tk.StringVar()
resultstyle_recon_var=Tk.StringVar()
Rcenter_recon_var=Tk.StringVar()
cbox_recon= ttk.Combobox()
label_recon=Tk.Label()
label_intg=Tk.Label()
Tk.Button(master=root,text='Quit',command=_quit,font=("Helvetica",20)).grid(row=figrowspan+1,column=0,rowspan=3,sticky=Tk.W,padx=10)

firstRowFrame = Tk.Frame(root)
firstRowFrame.grid(row=figrowspan+1,column=1,sticky=Tk.W)

Tk.Button(master=firstRowFrame,text='Please Choose Parameter File for Integration',command=folderSelector_parameter_intg,font=("Helvetica",12)).grid(row=25,column=2,pady=(20,0))
filePath_param_intg_en=Tk.Entry(master=firstRowFrame,width=100)
filePath_param_intg_en.grid(row=26,column=2,sticky=Tk.W)
filePath_param_intg_en.delete(0,"end")
filePath_param_intg_en.insert(0,"parameter file for integration")

Tk.Button(master=firstRowFrame, text='Open/Edit/Save Parameter File', command=open_parameterfile_intg,font=("Helvetica",12)).grid(row=26,column=3,padx=(10,10))

Tk.Button(master=firstRowFrame,text='Preview Data',command=linkToff_aysm,font=("Helvetica",14)).grid(row=29,column=2,pady=(10,10))


label_intg=Tk.Label(master=firstRowFrame, text=progress_integration, font=("Helvetica",12))
label_intg.grid(row=30,column=2,padx=10)
Tk.Button(master=firstRowFrame,text='Do Integration',font=("Helvetica",20),command=show_integration_progress).grid(row=34,column=2,rowspan=3,padx=10)
Tk.Button(master=firstRowFrame,text='Stop Integration',font=("Helvetica",20),command=stop_integration).grid(row=34,column=3,rowspan=3,padx=(10,100))


Tk.Button(master=firstRowFrame,text='Please Choose Parameter File for Checking Integration Results',command=folderSelector_parameter_intg_read,font=("Helvetica",12)).grid(row=38,column=2,pady=(30,10))
filePath_param_intg_read_en=Tk.Entry(master=firstRowFrame,width=100)
filePath_param_intg_read_en.grid(row=40,column=2,sticky=Tk.W)
filePath_param_intg_read_en.delete(0,"end")
filePath_param_intg_read_en.insert(0,"parameter file for checking integration")

Tk.Button(master=firstRowFrame,text='Please Press to Load Above Prameter File to Continue',command=readParams_intg,font=("Helvetica",12)).grid(row=43,column=2,pady=(10,10))

Tk.Button(master=firstRowFrame,text='Check Integration',command=createNewWindow_integration,font=("Helvetica",14)).grid(row=45,column=2,rowspan=3,pady=(10,20))

Tk.Button(master=firstRowFrame,text='Please Choose Parameter File for Reconstruction',command=folderSelector_parameter_recon,font=("Helvetica",12)).grid(row=48,column=2,pady=(30,10))
filePath_param_recon_en=Tk.Entry(master=firstRowFrame,width=100)
filePath_param_recon_en.grid(row=50,column=2,sticky=Tk.W)
filePath_param_recon_en.delete(0,"end")
filePath_param_recon_en.insert(0,"parameter file for reconstruction")

Tk.Button(master=firstRowFrame, text='Open/Edit/Save Parameter File', command=open_parameterfile_recon,font=("Helvetica",12)).grid(row=50,column=3,padx=(10,10))

label_recon=Tk.Label(master=firstRowFrame, text=progress_Reconstruction, font=("Helvetica",12))
label_recon.grid(row=51,column=2,padx=10)
Tk.Button(master=firstRowFrame,text='Run DT',font=("Helvetica",20),command=show_recon_progress).grid(row=52,column=2,rowspan=3,padx=10)
Tk.Button(master=firstRowFrame,text='Stop DT',font=("Helvetica",20),command=stop_progress).grid(row=52,column=3,rowspan=3,padx=(10,100))

Tk.Button(master=firstRowFrame,text='Please Choose Parameter File for Checking Reconstruction Results',command=folderSelector_parameter_wf,font=("Helvetica",12)).grid(row=55,column=2,pady=(30,10))
filePath_param_en=Tk.Entry(master=firstRowFrame,width=100)
filePath_param_en.grid(row=56,column=2,sticky=Tk.W)
filePath_param_en.delete(0,"end")
filePath_param_en.insert(0,"parameter file for checking results")

Tk.Button(master=firstRowFrame,text='Please Press to Load Above Parameter File to Continue',command=readParams,font=("Helvetica",12)).grid(row=58,column=2,pady=(10,10))


Tk.Button(master=firstRowFrame,text='Check LineOut',command=createNewWindow_lineout,font=("Helvetica",14)).grid(row=60,column=2,rowspan=3,pady=(20,10))
Tk.Button(master=firstRowFrame,text='Check Sinogram',command=createNewWindow_sino,font=("Helvetica",14)).grid(row=63,column=2,rowspan=3,pady=(10,10))
Tk.Button(master=firstRowFrame,text='Check Tomo',command=createNewWindow_tomo,font=("Helvetica",14)).grid(row=66,column=2,rowspan=3,pady=(10,10))
Tk.Button(master=firstRowFrame,text='Check Reconstruction',command=createNewWindow_recon,font=("Helvetica",14)).grid(row=69,column=2,rowspan=3,pady=(10,10))


Tk.mainloop()
