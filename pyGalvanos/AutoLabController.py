
import os
import clr
import sys
import numpy as np
import matplotlib.pyplot as plt
import time

hdwPath=R'C:\Program Files\Metrohm Autolab\Autolab SDK 2.1\Hardware Setup Files\uAutolabIII\HardwareSetup.xml'
sdkPath=R"C:\Program Files\Metrohm Autolab\Autolab SDK 2.1\EcoChemie.Autolab.Sdk"
adxPath=R"C:\Program Files\Metrohm Autolab\Nova 2.1\config\Adk.x"



print(os.path.dirname(sdkPath))
sys.path.append(os.path.dirname(sdkPath))
 
clr.AddReference("EcoChemie.Autolab.Sdk")
from EcoChemie.Autolab import Sdk as sdk
autolab = sdk.Instrument()
autolab.HardwareSetupFile = hdwPath
autolab.AutolabConnection.EmbeddedExeFileToStart=adxPath
autolab.Connect()

def CellOnOff(onoff):
        if onoff  == True:
            autolab.Ei.CellOnOff =sdk.EI.EICellOnOff.On
        elif onoff ==False:
            autolab.Ei.CellOnOff =sdk.EI.EICellOnOff.Off

def SetRefPotential(potential):
        autolab.Ei.Setpoint=potential*-1
        
SetRefPotential(0)  


def runCV_autolab( startPotential_mV, upperPotenial_mV, lowerPotential_mV, stopPotential_mV, cycles, scanrate_mV_s, plot=False):
        myProcedure = autolab.LoadProcedure(r"C:\Data\StandardCV.nox")
        myCommand = myProcedure.Commands["CV staircase"]
        myCommand.CommandParameters["Start potential (V)"].Value=startPotential_mV/1000
        myCommand.CommandParameters["Upper vertex potential (V)"].Value=upperPotenial_mV/1000
        myCommand.CommandParameters["Lower vertex potential (V)"].Value=lowerPotential_mV/1000
        #myCommand.CommandParameters["Step potential (V)"].Value=
        myCommand.CommandParameters["Number of stop crossings"].Value=cycles
        myCommand.CommandParameters["Stop potential (V)"].Value=stopPotential_mV/1000
        myCommand.CommandParameters["Scan rate (V/s)"].Value=scanrate_mV_s/1000
        
        
        Vs = [] 
        Cs = []
        lPot=0
        while (myProcedure.IsMeasuring):
            pot,cur=float(autolab.Ei.get_Potential()),float(autolab.Ei.Current)
            if lPot!=pot:
                Vs.append( pot )
                Cs.append( cur )
        if plot:
            plt.title('Autolab')
            plt.plot(Vs,Cs )
            plt.ylabel('Current (nA)')
            plt.xlabel('Bias (V vs Ag/AgCl) ')
            plt.show()
        return   np.array(Vs),np.array(Cs)
        
def runStrip_autolab(startVoltage_V,maxVoltage_V,slew_mV_s,  plot_Autolab=False):
    slew_V_s=slew_mV_s/1000.0

    samplesPerSec=200

    segmentTime1 = np.abs(maxVoltage_V-startVoltage_V)/slew_V_s
    segmentPoints =int( segmentTime1*samplesPerSec)
    S1=np.linspace(startVoltage_V,maxVoltage_V,segmentPoints)

    segmentTime2 = np.abs(maxVoltage_V-startVoltage_V)/slew_V_s
    segmentPoints =int( segmentTime2*samplesPerSec)
    S2=np.linspace(maxVoltage_V, startVoltage_V,segmentPoints)

    biasi = np.concatenate( [S1,S2])*-1
    
    autolab.Ei.CurrentRange =sdk.EI.EICurrentRange.CR09_10mA

    currents =[]
    rpotentials =[]
     
    SetRefPotential(0)
    CellOnOff(True)     
    time.sleep(5)
    for bias in biasi:
        SetRefPotential(bias)
        time.sleep(1.0/samplesPerSec)
        pot,cur=float(autolab.Ei.get_Potential()),float(autolab.Ei.Current)
        rpotentials.append(pot)
        currents.append(cur)

    SetRefPotential(0)
    CellOnOff(False)     

    if plot_Autolab:
        plt.title('Device')
        plt.ylabel('Current (nA)')
        plt.xlabel('Bias (V vs Ag/AgCl) ')
        plt.plot(-1*biasi,currents,label='autolab' )
        plt.legend()
            
    
    return biasi, currents,np.array(rpotentials) 