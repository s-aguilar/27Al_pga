import numpy as np
import pandas as pd

# Read in the data into dataframe
# df1 = pd.read_csv('Yields/P1/_P1.csv')
# df2 = pd.read_csv('Yields/P2/_P2.csv')
df3 = pd.read_csv('Yields/A1/_A1.csv')
df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet3')

pRun = df['Run'].values
pEproton = df['Ep (keV)'].values
pEproton = pEproton.round(2) # round to 2 decimal point


# Create dictionary, associating Run number to Proton energy
runToEproton = {pRun[i]:pEproton[i] for i in range(len(pRun))}

# dict_channels = {'p1':df1,'p2':df2,'a1':df3}
dict_channels = {'a1':df3}
# channels = ['p1','p2']#,'a1']
channels = ['a1']

for ch in channels:
    print('\n\nWORKING ON %s CHANNEL:'%ch)
    chan = dict_channels[ch]
    Eproton = []
    RunNum = []
    for ind, val in chan['Run'].iteritems():
        val = np.int64(val[4:])
        Eproton.append(runToEproton[val])
        RunNum.append(val)

    DetNum = []
    for ind, val in chan['Detector'].iteritems():
        val = np.int64(val[3:])
        DetNum.append(val)

    # Columns are not useful as type string, get rid of them and reindex them
    chan = chan.drop(columns = ['Run','Detector'])
    chan = chan.assign(Ep=pd.Series(Eproton,index=chan.index).values)
    chan = chan.assign(Run=pd.Series(RunNum,index=chan.index).values)
    chan = chan.assign(Detector=pd.Series(DetNum,index=chan.index).values)


    # Assign column as index
    # chan = chan.set_index(['Run','Ep'])
    chan = chan.set_index(['Run'])


    # Pre-Cleaning
    ## Drop data where the conditions are not met

    # Drop points with bad fits
    print('\n\nChecking if isValid is ever bad:')
    print(chan.query('IsValid != 1'))
    print('\nChecking if Status is ever bad:')
    print(chan.query('Status != 0'))

    chan.query('IsValid == 1',inplace=True)
    chan.query('Status == 0',inplace=True)


    # Drop points where error in yield is greater than the yield
    print('\nChecking if error in yield is greater than the yield:')
    print(chan.query('Yield < `Yield err`'))

    chan.query('Yield > `Yield err`',inplace=True)


    # Drop points where there are no counts
    print('\nChecking for points that have no counts:')
    print(chan.query('Area < 1'))

    chan.query('Area > 0',inplace=True)


    # Drop points that are from BG run if they snuck in (forgot to skip analysis of those runs)
    print('\nChecking for BG runs:')
    print(chan.query('Ep == 0'))

    chan.query('Ep != 0',inplace=True)


    # Pre-Cleaning complete, drop columns (no longer needed)
    chan.drop(columns = ['IsValid','Status','Q_int'], inplace=True)
    print(chan.head(5))


    # Save
    chan.to_csv('Yields/%s/%sYields.csv'%(ch.upper(),ch), sep=',')
    chan.to_excel('Yields/%s/%sYields.xlsx'%(ch.upper(),ch))

    # Make different sheets for each detector (quicker to plot and identify problems)
    with pd.ExcelWriter('Yields/%s/%sYieldsPerDetectors.xlsx'%(ch.upper(),ch)) as writer:
        for det in range(13):
            cut = (chan['Detector'] == det)
            tempdf = chan[cut]
            tempdf.to_excel(writer,sheet_name='%s'%det)

print('DONE!')




    #
    # # Read in the data into dataframe
    # # df1 = pd.read_csv('Yields/P1/_P1.csv')
    # # df2 = pd.read_csv('Yields/P2/_P2.csv')
    # df3 = pd.read_csv('Yields/A1/_A1.csv')
    # df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet3')
    #
    # pRun = df['Run'].values
    # pEproton = df['Ep (keV)'].values
    # pEproton = pEproton.round(2) # round to 2 decimal point
    #
    #
    # # Create dictionary, associating Run number to Proton energy
    # runToEproton = {pRun[i]:pEproton[i] for i in range(len(pRun))}
    #
    # Eproton = []
    # RunNum = []
    # for ind, val in df3['Run'].iteritems(): # Only for df3
    #     val = np.int64(val[4:])
    #     Eproton.append(runToEproton[val])
    #     RunNum.append(val)
    #
    # DetNum = []
    # for ind, val in df3['Detector'].iteritems(): # Only for df3
    #     val = np.int64(val[3:])
    #     DetNum.append(val)
    #
    # # Columns are not useful as type string
    # df3 = df3.drop(columns = ['Run','Detector'])
    #
    # # Append new columns (E proton) to dataframe, preserving the index
    # # df1 = df1.assign(Ep=pd.Series(Eproton,index=df1.index).values)
    # # df1 = df1.assign(RunNum=pd.Series(RunNum,index=df1.index).values)
    # # df2 = df2.assign(Ep=pd.Series(Eproton,index=df2.index).values)
    # # df2 = df2.assign(RunNum=pd.Series(RunNum,index=df2.index).values)
    # df3 = df3.assign(Ep=pd.Series(Eproton,index=df3.index).values)
    # df3 = df3.assign(Run=pd.Series(RunNum,index=df3.index).values)
    # df3 = df3.assign(Detector=pd.Series(DetNum,index=df3.index).values)
    #
    #
    # df3 = df3.set_index(['Run','Ep'])
    # # print(df3.head(20))
    # # exit()
    #
    # # df1.to_csv('Yields/P1/p1Yields.csv', sep=',')
    # # df1.to_excel('Yields/P1/p1Yields.xlsx')
    # #
    # # df2.to_csv('Yields/P2/p2Yields.csv', sep=',')
    # # df2.to_excel('Yields/P2/p2Yields.xlsx')
    #
    # df3.to_csv('Yields/A1/a1Yields.csv', sep=',')
    # df3.to_excel('Yields/A1/a1Yields.xlsx')
    #
    #
    # # Make different sheets for each detector (easier to plot and identify problems)
    # with pd.ExcelWriter('Yields/A1/a1YieldsPerDetectors.xlsx') as writer:
    #     for det in range(13):
    #         cut = (df3['Detector'] == det)
    #         tempdf = df3[cut]
    #         tempdf.to_excel(writer,sheet_name='%s'%det)
    # print('DONE!')
