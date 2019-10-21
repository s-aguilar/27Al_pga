import numpy as np
import pandas as pd


# Read in the data into dataframe
df1 = pd.read_csv('Yields/P1/_P1.csv')
df2 = pd.read_csv('Yields/P2/_P2.csv')
df3 = pd.read_csv('Yields/A1/_A1.csv')
df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet3')

pRun = df['Run'].values
pEproton = df['Ep (keV)'].values
pEproton = pEproton.round(2) # round to 2 decimal point


# Create dictionary, associating Run number to Proton energy
runToEproton = {pRun[i]:pEproton[i] for i in range(len(pRun))}

Eproton = []
RunNum = []
for ind, val in df1['Run'].iteritems():
    val = np.int64(val[4:])
    Eproton.append(runToEproton[val])
    RunNum.append(val)


# Append new columns (E proton) to dataframe, preserving the index
df1 = df1.assign(Ep=pd.Series(Eproton,index=df1.index).values)
df1 = df1.assign(RunNum=pd.Series(RunNum,index=df1.index).values)
df2 = df2.assign(Ep=pd.Series(Eproton,index=df2.index).values)
df2 = df2.assign(RunNum=pd.Series(RunNum,index=df2.index).values)
df3 = df3.assign(Ep=pd.Series(Eproton,index=df3.index).values)
df3 = df3.assign(RunNum=pd.Series(RunNum,index=df3.index).values)

df1.to_csv('Yields/P1/p1Yields.csv', sep=',')
df1.to_excel('Yields/P1/p1Yields.xlsx')

df2.to_csv('Yields/P2/p2Yields.csv', sep=',')
df2.to_excel('Yields/P2/p2Yields.xlsx')

df3.to_csv('Yields/A1/a1Yields.csv', sep=',')
df3.to_excel('Yields/A1/a1Yields.xlsx')
print('DONE!')
