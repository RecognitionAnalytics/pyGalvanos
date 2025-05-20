import numpy as np
import matplotlib.pyplot as plt 
import os
import pandas as pd
import glob
import json
from AutoEC_Lib import FindMolecules, DuckParams ,loadPine, loadGamry

#########################################Enter the folder with the CVs##############################################

allFiles=r'C:\Users\bashc\Downloads\scanrate\scan rate'

####################################################################################################################
outFolder= allFiles
allFiles = glob.glob(allFiles + "\\*.dta")
for file in allFiles:print(f"'{file}'")
 
 
 
CVs={}
cc=0
for dataDir in allFiles:
    #if the extension is .csv then load through a pine
    if os.path.splitext(dataDir)[1]=='.csv' or os.path.splitext(dataDir)[1]=='':
        cv=loadPine(dataDir)
    elif os.path.splitext(dataDir)[1].lower()=='.dta':
        cv=loadGamry(dataDir)
    plt.plot(cv.V,cv.C*1e6,label=os.path.basename(dataDir))
    CVs[os.path.basename(dataDir)]=cv

 
plt.title="All CVs"
plt.xlabel("E/V vs Ag/AgCl")
plt.ylabel("C (uA)") 
 
plt.show()


CVs={}
cc=0
tables ={}
for dataDir in allFiles:
    #if the extension is .csv then load through a pine
    if os.path.splitext(dataDir)[1]=='.csv' or os.path.splitext(dataDir)[1]=='':
        cv=loadPine(dataDir)
    elif os.path.splitext(dataDir)[1].lower()=='.dta':
        cv=loadGamry(dataDir)
     
    filename = os.path.basename(dataDir)
    params=DuckParams(filename,cv,RedoxConcentration_M=5e-3)
    plt.savefig(outFolder + "\\" + filename + '.jpg')
    plt.show() 
    tables[filename]=params
 
 
t=pd.DataFrame(tables).transpose()
t.to_csv(outFolder + "\\tables.csv") 