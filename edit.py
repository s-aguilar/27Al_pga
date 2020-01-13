'''
This script reads in the yield files produced from `xxyields_pa.C` as well as the
experimental online logbook file. A pandas DataFrame is used, and modified to
have the yields, the efficiency corrected yields, run number, and beam energy.
The DataFrame is then `cleaned` for bad data and when finalized it is then
saved to an excel file
'''
import numpy as np
import pandas as pd

# Read online analysis logbook
df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet3')


# Read in efficiencies
dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')


pRun = df['Run'].values
pEproton = df['Ep (keV)'].values
pEproton = pEproton.round(0) # round to 2 decimal points


# Create dictionary, associating Run number to Proton energy
runToEproton = {pRun[i]:pEproton[i] for i in range(len(pRun))}


# channels = ['a1','p1','p2']
channels = ['a1']

for ch in channels:

    print('\n\nWORKING ON %s CHANNEL:'%ch)

    # Read in the data into dataframe
    dtypeDict = {'Yield':np.float,'Yield err':np.float,'Area':np.float,'Area err':np.float,
                'X2NDF':np.float,'isValid':np.int,'Status':np.int,'Q int':np.int}
    chan = pd.read_csv('Yields/%s/_%s.csv'%(ch.upper(),ch.upper()),dtype=dtypeDict)
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
    # chan = chan.set_index(['Run'])


    # Efficiency correct the yields and append to DataFrame
    eff = dfeff[ch].values
    angle = dfeff['Angle'].values
    chan = chan.assign(Angle=pd.Series(angle[0:len(chan['Yield'].values)],index=chan.index).values)


    # Convert from pulses to Coulombs and apply efficieny corrections
    q_e = 1.602e-19
    scale = 1e-8    # 10^-8 C/pulse
    q_corr = scale/(q_e)
    yield_effcor = chan['Yield'].values/ q_corr / eff[0:len(chan['Yield'].values)]
    yield_err_effcor = chan['Yield err'].values/ q_corr / eff[0:len(chan['Yield'].values)]


    # Update the DataFrame with efficiency corrected yields
    chan = chan.assign(Yield_effcor=pd.Series(yield_effcor,index=chan.index).values)
    chan = chan.assign(Yield_err_effcor=pd.Series(yield_err_effcor,index=chan.index).values)
    chan = chan.rename(columns={'Yield_effcor': 'Yield effcor', 'Yield_err_effcor': 'Yield err effcor'})


    # Explicitly drop BG runs
    chan.query('Run != 863',inplace=True)   # Explicit BG run
    chan.query('Run != 939',inplace=True)   # Explicit BG run


    # There are runs with identical beam energies, identify the runs and
    # average the yields per detector, update the data frame
    energySet = set(chan['Ep'].values)
    listOfDuplicateNRGRuns = []
    for nrg in energySet:
        nrgMask = chan.query('Ep == %f'%nrg)
        if len(nrgMask) > 13:
            listOfDuplicateNRGRuns.append(set(nrgMask['Run']))
    print(listOfDuplicateNRGRuns)

    # Run numbers with duplicate energies have been identified, now average
    # per detector, and reassign them back into DataFrame with new run numbers
    # starting from 2000
    newRunNum = 2000
    for dupeRun1,dupeRun2 in listOfDuplicateNRGRuns:

        tempdf = chan.query('Run == %f or Run == %f'%(dupeRun1,dupeRun2))
        # print(tempdf)

        # Custom function to correctly compute the average yield (weighted by charge)
        avgY = lambda x : np.average(x, weights=tempdf.loc[x.index,'Q int'])

        # Custom function to correctly compute the average of error (square root
        # of the sum in quadrature)
        # avgYErr = lambda x : np.sqrt(sum(np.array(x)**2))
        avgYErr = lambda x : np.sqrt(np.average(x, weights=tempdf.loc[x.index,'Q int']))

        # Custom function to give a new run number starting from 2000
        runnum = lambda x : newRunNum

        # Dictionary containing aggregation functions
        tempdict = {'Yield':avgY,'Yield err':avgYErr,
        'Area':'mean','Area err':avgYErr,'sig1':'mean','X2NDF':'mean',
        'IsValid':'mean','Status':'mean','Q int':'mean','Ep':'first',
        'Run':runnum,'Angle':'first','Yield effcor':avgY,
        'Yield err effcor':avgYErr
        }

        # New DataFrame containing detector averaged quantities (yield, area, etc)
        # as well as old quantities
        tempdf = tempdf.groupby('Detector').agg(tempdict)
        tempdf.reset_index(inplace=True)
        # print(tempdf)

        # Combine the new averaged DataFrame with the old
        chan = pd.concat([chan,tempdf],sort=True, axis=0)

        # Drop the duplicate energy runs
        chan.query('Run != %f'%dupeRun1,inplace=True)
        chan.query('Run != %f'%dupeRun2,inplace=True)

        newRunNum+=1


    # Pre-Cleaning (Drop data where these conditions are not met)

    # Drop points with bad fits
    # IsValid == 1 -> Good Fit
    # IsValid == 0 -> Bad Fit
    print('\n\nChecking if fit `isValid` is ever bad:')
    print(chan.query('IsValid != 1'))
    print('\n\nChecking if fit `Status` is ever bad:')
    print(chan.query('Status != 0'))

    chan.query('IsValid == 1',inplace=True)
    chan.query('Status == 0',inplace=True)


    # Drop points where error in yield is greater than the yield
    print('\n\nChecking if error in yield is greater than the yield:')
    print(chan.query('Yield < `Yield err`'))

    chan.query('Yield > `Yield err`',inplace=True)


    # Drop points where there are no counts
    print('\n\nChecking for points that have no counts:')
    print(chan.query('Area < 1'))

    chan.query('Area > 0',inplace=True)


    # Drop points that are from a BG run if they snuck in (were not explicitly
    # dropped)
    print('\n\nChecking for BG runs that may have snuck in:')
    print(chan.query('Ep == 0'))
    chan.query('Ep != 0',inplace=True)


    # Example how to drop all data of a specific Run number
    # chan.query('Run != 9999',inplace=True)


    # Average the yield for each detector if there are multiple runs with the
    # same energy.


    # Pre-Cleaning complete, drop columns (no longer needed)
    chan.drop(columns = ['IsValid','Status'], inplace=True)


    print('\n\nFinal form of DataFrame:')
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
