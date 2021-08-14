# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 17:03:49 2021

@author: GrossAR

This program reads in csv tables listing karyotype results and separates the normal
and abnormal lines.
"""


##############################################################################
### 1. Import Libraries

import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
np.set_printoptions(edgeitems=3, infstr='inf',linewidth=200, nanstr='nan', precision=8,suppress=False, threshold=1000, formatter=None)
#import scipy
##############################################################################
### 2. Define functions



##############################################################################
### 3. Import csv files

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Karyotype_lists')
filesToImport = os.listdir()

currentDF = pd.read_csv(filesToImport[0])
currentDF = currentDF.iloc[0:20,0:4]

##############################################################################
### 4. Formatting

###### 4.1 - Remove data bars

len(currentDF)
currentDF.loc('Passage')




cellLine = currentDF['Cell Line']

dateRow = []
for element in cellLine:
    #print(element)
    dateRow.append('Date' in element)

pd.concat([currentDF, pd.DataFrame(dateRow)], axis = 1, ignore_index=True)

    ret_value.append()
ret_value = 'Date' in cellLine
.head()