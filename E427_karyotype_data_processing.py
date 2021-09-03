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
from datetime import datetime
import math

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 1000)
np.set_printoptions(edgeitems=3, infstr='inf',linewidth=200, nanstr='nan', precision=8,suppress=False, threshold=1000, formatter=None)
pd.options.mode.chained_assignment = None  # default='warn'
#import scipy
##############################################################################
### 2. Define functions
###### 2.1 - Converts date elements into a standard date format
def convertDateFormat(dataframe, columnName, returnErrors = False):
    newKaryoDateCol = []
    counter = 0
    rowError =[]
    for date in dataframe[columnName]:
        try:
            newKaryoDateCol.append(datetime.strptime(date.replace('/','-'), '%M-%d-%Y').date())
        except:
                try:
                    newKaryoDateCol.append(datetime.strptime(date.replace('/','-'), '%M-%d-%y').date())
                except:
                        print('" ' + date + ' " in row ' + str(counter) + ' failed to format')
                        rowError.append(counter)
                        newKaryoDateCol.append(datetime(1, 1, 1).date())
              
        counter += 1

    if returnErrors == False:
        return(newKaryoDateCol)
    else:
        return(rowError)
    
def grep(pattern, aList):
    aList = list(aList)
    start = 0
    output = []
    while True:
        try:
            pos = aList.index(pattern, start)
            output.append(pos)
            start = pos+1
        except:
            return(output)
            break
    
###### 2.2 - Bins ages
def binAges(dataframe):  ### Returns a new column as a new dataframe
    ageBin = pd.DataFrame(np.arange(len(dataframe))*0)
    ageBin[dataframe['Age'].between(0   , 20, inclusive = True).tolist()] = 1
    ageBin[dataframe['Age'].between(20.1, 40, inclusive = True).tolist()] = 2
    ageBin[dataframe['Age'].between(40.1, 60, inclusive = True).tolist()] = 3
    ageBin[dataframe['Age'].between(60.1, 80, inclusive = True).tolist()] = 4
    ageBin[dataframe['Age'].between(80.1, 200, inclusive = True).tolist()] = 5
    ageBin = ageBin.loc[:,0].tolist()
    return(ageBin)

def binAges2(dataframe, binSize):  ### Returns a new column as a new dataframe
    numberOfBins = math.ceil(100/ binSize)
    ageBin = pd.DataFrame(np.arange(len(dataframe))*0)
    for bin in range(0,numberOfBins):
        ageBin[dataframe['Age'].between(bin*binSize+0.01   , (bin+1)*binSize, inclusive = True).tolist()] = bin+1
    print(pd.DataFrame(ageBin).value_counts().sort_index())
    ageBin = ageBin.loc[:,0].tolist()

    return(ageBin)

##############################################################################
### 3. Import csv files

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Karyotype_lists')
filesToImport = os.listdir()
filesToImport2 = []
for file in filesToImport:
    print(file)
    if '.csv' in file:
        filesToImport2.append(file)
    else:
        pass

filesToImport = filesToImport2
del(filesToImport2)

currentDF = pd.read_csv(filesToImport[0]).iloc[0:20,0:5]

##############################################################################
### 4. Formatting

###### 4.1 - Join each csv into a shared DF
completeDF = pd.read_csv(filesToImport[0]).iloc[0:1,0:5]
#completeDF = completeDF.iloc[0]

for file in filesToImport:
    print(file)
    newDF = pd.read_csv(file).iloc[:,0:5]
    print(newDF.shape)
    completeDF = completeDF.append(newDF)

len(completeDF)
completeDF = completeDF.dropna()

###### 4.2 - Generate new column from date received

dataframe = completeDF[0:300]
repeatCounter = 0
newCol = []
for rowNum in range(0,len(dataframe)):
    rowCurrent = dataframe.iloc[rowNum]
    try:
        if 'Date received' in rowCurrent[0]:
            latestDate = rowCurrent[0].split()[2].replace('/','-')
            latestDate = datetime.strptime(latestDate,'%M-%d-%y').date()
        else:
            pass
    except:
        print('Error in row ' + str(rowNum) + ': \n' + rowCurrent)
    newCol.append(latestDate)
    #print(rowCurrent)

dataframe['Date_Recieved'] = newCol

###### 4.3 - Generate new column from date received

dataframe = dataframe[dataframe['Passage'].notnull()]

###### 4.4 - Convert Date of Karyotype to standard date format
'''
newKaryoDateCol = []
counter = 0
rowError =[]
for date in dataframe['Date of Karyotype ']:
    try:
        newKaryoDateCol.append(datetime.strptime(date.replace('/','-'), '%M-%d-%Y').date())
    except:
        try:
            newKaryoDateCol.append(datetime.strptime(date.replace('/','-'), '%M-%d-%y').date())
        except:
            print(date + ' in row ' + str(counter) + ' failed to format')
            rowError.append(counter)
    counter += 1
'''
newKaryoDateCol = convertDateFormat(dataframe, 'Date of Karyotype ')
dataframe['Date of Karyotype '] = newKaryoDateCol

###### 4.5 - Add a column marking normal or abnormal
newColumn = []
for rowNum in range(0,len(dataframe)):
    rowCurrent = dataframe.iloc[rowNum]
    if 'Abnormal' in rowCurrent['Results'] :
        newColumn.append('Abnormal')
    else:
        newColumn.append('Nominal')

dataframe['Normality'] = newColumn


###### 4.6 - Add a column marking expanded or non

print(dataframe['Parent cell type'].unique())
counter = 0
newColumn = []

for rowNum in range(0,len(dataframe)):
    rowCurrent = dataframe.iloc[rowNum]
    if 'Fibroblast' in rowCurrent['Parent cell type'] :
        newColumn.append('Expanded')
    
    elif 'LCL' in rowCurrent['Parent cell type'] :
        newColumn.append('Expanded')
    
    elif 'CJE' in rowCurrent['Parent cell type'] :
        newColumn.append('Expanded')

    elif 'PBMC' in rowCurrent['Parent cell type'] :
        newColumn.append('Unex')
        
    else:
        print('Error: unrecognized cell type in line ' + str(counter) + ': ' + rowCurrent['Parent cell type'])

    counter += 1


dataframe['Expansion'] = newColumn

###### 4.7 - Add a parent line name
dataframe['ParentName'] = dataframe['Cell Line'].str.split('-', n=1, expand=True).iloc[:,0]
cols = dataframe.columns.tolist()
cols = cols[0:1] + ['ParentName'] + cols[1:-1]
dataframe = dataframe[cols]

len(dataframe['ParentName'].unique())

dataframe.iloc[50:60]

##############################################################################
### 5. Save / export new dataframe
###### 5.1 - Add a column marking expanded or non

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates')
dataframe.to_csv('karyotypes-FULL.csv')
dataframe.to_excel('karyotypes-FULL.xls', index=False)


##############################################################################
### 6. Import additional metadata
###### 6.1 - Sort through and isolate csv files for processing
os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates/Karyotype_lists_by_line')
filesToImport = os.listdir()
filesToImport2 = []
for file in filesToImport:
    print(file)
    if '.csv' in file:
        filesToImport2.append(file)
    else:
        pass

filesToImport = filesToImport2
del(filesToImport2)


###### 6.2 - Join each csv into a shared DF

columnNames = ['Donor #', 'Cell Line', 'Type', 'Passage-reprogramming', 'Age', 'Sex', 'Source', 'Passage-Karyo-1', 'Date-Karyo-1', 'Results-Karyo-1', 'Group-Karyo-1', 'Passage-Karyo-2', 'Date-Karyo-2', 'Results-Karyo-2', 'Group-Karyo-2', 'Passage-Karyo-3', 'Date-Karyo-3', 'Results-Karyo-3', 'Group-Karyo-3', 'Passage-Karyo-4', 'Date-Karyo-4', 'Results-Karyo-4', 'Group-Karyo-4', 'Passage-Karyo-5', 'Date-Karyo-5', 'Results-Karyo-5', 'Group-Karyo5']
metadata = pd.DataFrame(columns = columnNames)

for file in filesToImport:
    print(file)
    newDF = pd.read_csv(file)
    print(newDF.shape)
    metadata = metadata.append(newDF)

### For troubleshooting:
#metadata = pd.read_csv(filesToImport[4]).iloc[0:1,]

##############################################################################
### 7. Create a master list of every karyotype from the metadata
###### 7.1 - Loop throug and check how many karyotypes each row has, then give each a new line
### Takes around 1.2 s/100 rows
columnNames = ['Donor #', 'Cell Line', 'Type', 'Passage-reprogramming', 'Age', 'Sex', 'Source', 'Passage-Karyo', 'Date-Karyo', 'Result', 'Normality', 'KaryoNum']
newDF = pd.DataFrame(columns= columnNames)

for rowNum in range(0,len(metadata)):
    rowCurrent = metadata.iloc[rowNum]   # Isolate the current row to process 
    newBlock = rowCurrent[0:11]
    newBlock['KaryoNum'] = 1

    newBlock = pd.DataFrame([newBlock.tolist()], columns= columnNames)
    if math.isnan(rowCurrent['Passage-Karyo-2'])==False:
        rowNew = rowCurrent[0:7]
        rowNew = rowNew.append(rowCurrent[11:15])
        rowNew['KaryoNum'] = 2
        rowNew = pd.DataFrame([rowNew.tolist()], columns= columnNames)
        newBlock = newBlock.append(rowNew)
    if math.isnan(rowCurrent['Passage-Karyo-3'])==False:
        rowNew = rowCurrent[0:7]
        rowNew = rowNew.append(rowCurrent[15:19])
        rowNew['KaryoNum'] = 3
        rowNew = pd.DataFrame([rowNew.tolist()], columns= columnNames)
        newBlock = newBlock.append(rowNew)
    if math.isnan(rowCurrent['Passage-Karyo-4'])==False:
        rowNew = rowCurrent[0:7]
        rowNew = rowNew.append(rowCurrent[19:23])
        rowNew['KaryoNum'] = 4
        rowNew = pd.DataFrame([rowNew.tolist()], columns= columnNames)
        newBlock = newBlock.append(rowNew)
    if math.isnan(rowCurrent['Passage-Karyo-5'])==False:
        rowNew = rowCurrent[0:7]
        rowNew = rowNew.append(rowCurrent[23:27])
        rowNew['KaryoNum'] = 5
        rowNew = pd.DataFrame([rowNew.tolist()], columns= columnNames)
        newBlock = newBlock.append(rowNew)
    newDF = newDF.append(newBlock)

karyoFull = newDF

###### 7.2.1 - Remove Source column
colToKeep = karyoFull.columns.tolist()
colToKeep
del colToKeep[6]
karyoFull = karyoFull[colToKeep]

###### 7.3 - Convert Date of Karyotype to standard date format
karyoFull['Date-Karyo'] = (karyoFull['Date-Karyo']).apply(str)
newDates = convertDateFormat(karyoFull, 'Date-Karyo')
karyoFull['Date-Karyo'] = newDates

###### 7.4 - Add a column marking expanded or non
print(karyoFull['Type'].unique())
counter = 0
newColumn = []
for rowNum in range(0,len(karyoFull)):
    rowCurrent = karyoFull.iloc[rowNum]
    if 'Fibroblast' in rowCurrent['Type'] :
        newColumn.append('Expanded')
    
    elif 'LCL' in rowCurrent['Type'] :
        newColumn.append('Expanded')
    
    elif 'Epithelial' in rowCurrent['Type'] :
        newColumn.append('Expanded')
    
    elif 'Adipose' in rowCurrent['Type'] :
        newColumn.append('Expanded')

    elif 'PBMC' in rowCurrent['Type'] :
        newColumn.append('Unexpanded')
        
    else:
        print('Error: unrecognized cell type in line ' + str(counter) + ': ' + rowCurrent['Parent cell type'])

    counter += 1


karyoFull['Expansion'] = newColumn


###### 7.5 - Add a column with the parent cell line name
parentNames = []
for rowNum in range(0,len(karyoFull)):
    parentName = karyoFull['Cell Line'].iloc[rowNum].split('-')[0]
    parentNames.append(parentName)

karyoFull['Parent Line'] = parentNames

###### 7.6 - Fill in donor number

uniqueParentNames = list(dict.fromkeys(parentNames))

donorNum = []
for rowNum in range(0,len(karyoFull)):
    parentName = parentNames[rowNum]
    newDonorNum = uniqueParentNames.index(parentName) +1
    donorNum.append(newDonorNum)

karyoFull['Donor #'] = donorNum


###### 7.8 - Add a row index
rowIndex = list(range(1,len(karyoFull)+1))
karyoFull['RowNum'] = rowIndex

###### 7.7 - Reorder Columns
newOrder = ['RowNum','Donor #',
 'Parent Line', 
 'Cell Line',
 'Type',
 'Expansion',
 'Age',
 'Sex',
 'Passage-reprogramming',
 'Passage-Karyo',
 'Date-Karyo',
 'KaryoNum',
 'Result',
 'Normality']

karyoFull = karyoFull[newOrder]

#colNames = metadata.columns.tolist()
#colKeep = colNames[1:3] + [colNames[4]] + colNames[7:10] + colNames[11:14] + colNames[15:18] + colNames[19:22] + colNames[23:26]
#metadata = metadata[colKeep]

###### 8.0 - Export

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates')
karyoFull.to_excel('karyotypes-FULL.xls', index=False)

###### 9.0 Reimport if necessary
karyoImport = pd.read_excel('karyotypes-FULL.xls')


###### 10.0 - Extract key stats
### 10.1 - subset data by source cell type

epi = karyoImport[karyoImport['Type']=='Epithelial']
adi = karyoImport[karyoImport['Type']=='Adipose']
fib = karyoImport[karyoImport['Type']=='Fibroblast']
mfib = karyoImport[karyoImport['Type']=='Mod-Fibroblast']
lcl = karyoImport[karyoImport['Type']=='LCL']
pbmc =karyoImport[karyoImport['Type']=='PBMC']
mpbmc = karyoImport[karyoImport['Type']=='Mod-PBMC']

### 10.2 - count parent lines

len(list(dict.fromkeys(epi['Parent Line'])))
len(list(dict.fromkeys(adi['Parent Line'])))
len(list(dict.fromkeys(fib['Parent Line'])))
len(list(dict.fromkeys(mfib['Parent Line'])))
len(list(dict.fromkeys(lcl['Parent Line'])))
len(list(dict.fromkeys(pbmc['Parent Line'])))
len(list(dict.fromkeys(mpbmc['Parent Line'])))

### 10.3 - count clone lines

len(list(dict.fromkeys(epi['Cell Line'])))
len(list(dict.fromkeys(adi['Cell Line'])))
len(list(dict.fromkeys(fib['Cell Line'])))
len(list(dict.fromkeys(mfib['Cell Line'])))
len(list(dict.fromkeys(lcl['Cell Line'])))
len(list(dict.fromkeys(pbmc['Cell Line'])))
len(list(dict.fromkeys(mpbmc['Cell Line'])))

### 10.4 - count karyotypes

len(epi)
len(adi)
len(fib)
len(mfib)
len(lcl)
len(pbmc)
len(mpbmc)

### 10.1 - count expanded line stats

karyoEx = karyoImport[karyoImport['Expansion'] == 'Expanded']
len(karyoEx)
np.count_nonzero(karyoEx['Normality']=='Abnormal')

karyoUnex = karyoImport[karyoImport['Expansion'] == 'Unexpanded']
len(karyoUnex)
np.count_nonzero(karyoUnex['Normality']=='Abnormal')


len(list(dict.fromkeys(karyoUnex['Cell Line'])))


cellType = 'Fibroblast'
cellType = 'PBMC'
cellType = 'Epithelial'
(list(np.where(karyoImport['Type'] == cellType))[0].tolist())
len(list(np.where(newDF['Type'] == cellType))[0].tolist())

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates')
karyoExAb = karyoEx[karyoEx['Normality'] == 'Abnormal']
karyoExAb.to_excel('karyotypes-ExAb.xls', index=False)

karyoUnexAb = karyoUnex[karyoUnex['Normality'] == 'Abnormal']
karyoUnexAb.to_excel('karyotypes-UnexAb.xls', index=False)


'''
##############################################################################
### 8. Format additional metadata
###### 8.1 - Add a parent line name
metadata['ParentName'] = metadata['Cell Line'].str.split('-', n=1, expand=True).iloc[:,0]
cols = metadata.columns.tolist()
cols = cols[0:2] + ['ParentName'] + cols[2:-1]
metadata = metadata[cols]

len(metadata['ParentName'].unique())

len(metadata)
metadata = metadata[metadata['Donor #'].notnull()]
len(metadata)


###### 8.3 - Create or add to a list of dictionaries

ageDict = {}
for rowNum in range(0,len(metadata)):
    line = metadata.iloc[rowNum]
    parentName = line['ParentName']
    age = line['Age']
    ageDict[parentName] = age

len(ageDict)

###### 8.2 - Create or add to a list of dictionaries
counter = 0
for rowNum in range(0,len(metadata)):
    rowCurrent = metadata.iloc[rowNum]
    if 'CS' in rowCurrent['ParentName'] :
        counter +=1
    else:
        pass


for rowNum in range(0,len(metadata)):
    line = metadata.iloc[rowNum]
    parentName = line['ParentName']
    try:
        parentName = parentName.split('CS')[-1]
    except:
        pass
    
    age = line['Age']
    ageDict[parentName] = age


os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates')
metadata.to_csv('metadata-FULL.csv')

##############################################################################
### 9. Add metadata to main data frame
###### 9.1 - Add age entries

ageColumn = []
found = 0
notFound = []
for rowNum in range(0,len(dataframe)):
    parentName = dataframe.iloc[rowNum]['ParentName']
    try:
        age = ageDict[parentName]
        found += 1
    except:
        #print(parentName + ' not found')
        age = 'NaN'
        notFound.append(parentName)
    ageColumn.append(age)
    
found
len(notFound)

dataframe['Age'] = ageColumn

###### 9.1 - Add age entries



##############################################################################
### 10. Extract age and passage number stats by group
###### 10.0 - Troubleshooting and data validation tests

karyoFull['Type'].value_counts()

test = binAges2(karyoFull[karyoFull['Type']=='Fibroblast'],10)
test = binAges2(karyoFull[karyoFull['Type']=='PBMC'],10)
test = binAges2(karyoFull[karyoFull['Type']=='LCL'],10)
test = test[test['Age']<=10]
test = binAges2(karyoFull[karyoFull['Type']=='Epithelial'],10)
test = binAges2(karyoFull[karyoFull['Type']=='Mod-Fibroblast'],10)
test = binAges2(karyoFull[karyoFull['Type']=='Mod-PBMC'],10)


###### 10.1 - Separate out Expanded an unexpanded
dataEx = karyoFull[karyoFull['Expansion']=='Expanded']
dataUnex = karyoFull[karyoFull['Expansion']=='Unex']

abnormalLines = karyoFull[karyoFull['Normality']=='Abnormal']
abnormalLines['Type'].value_counts()

dataExAb = abnormalLines[abnormalLines['Expansion']=='Expanded']
dataUnexAb = abnormalLines[abnormalLines['Expansion']=='Unex']

#dataExAb.to_excel('ExAb.xls')
#dataUnexAb.to_excel('UnexAb.xls')

###### 10.2 - Assign age bins

dataEx['AgeBin'] = binAges(dataEx)
dataEx['AgeBin'].value_counts().sort_index()
dataEx['AgeBin'] = binAges2(dataEx,20)

dataUnex['AgeBin'] = binAges(dataUnex)
dataUnex['AgeBin'].value_counts().sort_index()
dataUnex['AgeBin'] = binAges2(dataUnex,20)

dataExAb['AgeBin'] = binAges(dataExAb)
dataExAb['AgeBin'].value_counts().sort_index()
dataExAb['AgeBin'] = binAges2(dataExAb,20)

dataUnexAb['AgeBin'] = binAges(dataUnexAb)
dataUnexAb['AgeBin'].value_counts().sort_index()
dataUnexAb['AgeBin'] = binAges2(dataUnexAb,20)

###### 10.3 - Separate out first karyotypes from repeats
ExK1 = dataEx[dataEx['KaryoNum']==1]
ExK2 = dataEx[dataEx['KaryoNum']>1]
UnexK1 = dataUnex[dataUnex['KaryoNum']==1]
UnexK2 = dataUnex[dataUnex['KaryoNum']>1]

len(ExK1)                           # Number of total expanded first karyotypes
len(ExK2)                           # Number of total expanded subsequent karyotypes
len(UnexK1)                         # Number of total unexpanded first karyotypes
len(UnexK2)                         # Number of total unexpanded subsequent karyotypes

ExK1['Normality'].describe()        # Number of abnormal expanded first karyotypes
ExK2['Normality'].describe()        # Number of abnormal expanded subsequent karyotypes
UnexK1['Normality'].describe()      # Number of abnormal unexpanded first karyotypes
UnexK2['Normality'].describe()      # Number of abnormal unexpanded subsequent karyotypes

###### 10.4 - Separate by passage number for first karyotypes
ExLowP   = ExK1[ExK1['Passage-Karyo']<9]
ExHighP  = ExK1[ExK1['Passage-Karyo']>=9]
UnexLowP = UnexK1[UnexK1['Passage-Karyo']<9]
UnexHighP= UnexK1[UnexK1['Passage-Karyo']>=9]

len(ExLowP)                           # Total expanded karyotypes from passage 1-8
len(ExHighP)                          # Total expanded karyotypes from passage 9+
len(UnexLowP)                         # Total unexpanded karyotypes from passage 1-8
len(UnexHighP)                        # Total unexpanded karyotypes from passage 9+

ExLowP['Normality'].describe()        # Abnorm expanded karyotypes from passage 1-8
ExHighP['Normality'].describe()       # Abnorm expanded karyotypes from passage 9+
UnexLowP['Normality'].describe()      # Abnorm unexpanded karyotypes from passage 1-8
UnexHighP['Normality'].describe()     # Abnorm unexpanded karyotypes from passage 9+

###### 10.4 - Separate by passage number for subsequent karyotypes
ExLowP   = ExK2[ExK2['Passage-Karyo']<15]
ExHighP  = ExK2[ExK2['Passage-Karyo']>=15]
UnexLowP = UnexK2[UnexK2['Passage-Karyo']<15]
UnexHighP= UnexK2[UnexK2['Passage-Karyo']>=15]

len(ExLowP)                           # Total expanded karyotypes from passage 1-14
len(ExHighP)                          # Total expanded karyotypes from passage 15+
len(UnexLowP)                         # Total unexpanded karyotypes from passage 1-14
len(UnexHighP)                        # Total unexpanded karyotypes from passage 15+

ExLowP['Normality'].describe()        # Abnorm expanded karyotypes from passage 1-14
ExHighP['Normality'].describe()       # Abnorm expanded karyotypes from passage 15+
UnexLowP['Normality'].describe()      # Abnorm unexpanded karyotypes from passage 1-14
UnexHighP['Normality'].describe()     # Abnorm unexpanded karyotypes from passage 15+






len(dataEx)

len(dataUnex)

os.chdir('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E427 - Ideogram updates')
karyoFull.to_excel('Karyotypes-FULL.xls')
'''