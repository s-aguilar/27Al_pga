'''
This script reads in the yield files produced from `xxyields_pa.C` as well as
the experimental online logbook file. A pandas DataFrame is used, and modified
to have the yields, the efficiency corrected yields, run number, and beam
energy. The DataFrame is then `cleaned` for bad data and when finalized it is
then saved to an excel file, for subsequent analysis and plotting in
crossSectionV2.py
'''
import time
import numpy as np
import pandas as pd


# VERBOSITY
#   0 = off
#   1 = some
#   2 = all
verbose = 0



# Read online analysis logbook
df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet3')


# Read in efficiencies
dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')


pRun = df['Run'].values
pEproton = df['Ep (keV)'].values
pEproton = pEproton.round(1) # round to 2 decimal points
# print(pEproton,len(pEproton))
# pEproton = pEproton.round(0) # round to 2 decimal points
# pEproton = pEproton.astype(np.int)
# print(pEproton,len(set(pEproton)))
# exit()


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


    # Overlapping scans in data seem to appear to have an energy spread, try to
    # correct for it manually by shifting entire scans by some arbitrary amount
    # found by eye

    # Scan 5 is 0.5 keV shifted above scan 4
    scan5 = ((chan['Run'] > 807) & (chan['Run'] < 953))  # run numbers with endpoints exclusive runs 808-952
    # Scan 6 is 0.5 keV shifted relative above scan 5
    scan6 = ((chan['Run'] > 959) & (chan['Run'] < 1132))  # run numbers with endpoints exclusive runs 960-1131
    # print(chan.loc[scan6, 'Ep'].head())
    chan.update(chan.loc[scan5, 'Ep'] - 0)
    chan.update(chan.loc[scan6, 'Ep'] - .5)
    # print(chan.loc[scan6, 'Ep'].head()
    # exit()

    # Example how to drop all data of a specific Run number
    # chan.query('Run != 9999',inplace=True)


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


    # Beginning DataFrame
    start_time = time.time()
    print('\n\nFormat of beginnning DataFrame:')
    print(chan.head())
    print('Size:',len(chan),'rows')


    if verbose > 1:
        # IsValid == 1 -> Good Fit
        # IsValid == 0 -> Bad Fit
        print('\n\nChecking if fit `isValid` is ever bad:')
        print(chan.query('IsValid != 1'))
        print('\n\nChecking if fit `Status` is ever bad:')
        print(chan.query('Status != 0'))
    # Drop points with bad fits
    chan.query('IsValid == 1',inplace=True)
    chan.query('Status == 0',inplace=True)


    if verbose > 0:
        print('\n\nChecking if error in yield is greater than the yield:')
        print(chan.query('Yield < `Yield err`'))
        print('Number of points lost: ',len(chan.query('Yield < `Yield err`')))
    # Drop points where error in yield is greater than the yield
    chan.query('Yield > `Yield err`',inplace=True)


    if verbose > 0:
        print('\n\nChecking for points with low statistics:')
        print(chan.query('Area < 70'))
        print('Number of points lost: ',len(chan.query('Area < 70')))
    # Drop points with low statistics
    chan.query('Area > 70',inplace =True)


    if verbose > 1:
        print('\n\nChecking for BG runs that may have snuck in:')
        print(chan.query('Ep == 0'))
    # Drop points that are from a BG run if they snuck in (were not explicitly
    # dropped)
    chan.query('Ep != 0',inplace=True)


    # There are runs with identical beam energies, identify the runs and
    # average the yields per detector, update the data frame
    energySet = set(chan['Ep'].values)
    listOfDuplicateNRGRuns = []
    for nrg in energySet:
        nrgMask = chan.query('Ep == %f'%nrg)
        if len(nrgMask) > 13:
            listOfDuplicateNRGRuns.append(set(nrgMask['Run']))
    # print(listOfDuplicateNRGRuns)

    # Run numbers with duplicate energies have been identified, now average
    # per detector, and reassign them back into DataFrame with new run numbers
    # starting from 2000
    newRunNum = 2000
    for dupeRunSet in listOfDuplicateNRGRuns:

        # Unpack the set into a list to access the elements
        dupeRunList = [x for x in dupeRunSet]

        customQuery = ''
        for runs in dupeRunList:
            customQuery += 'Run == %d or '%runs

        # Remove extraneous ' or ' for last appending query
        customQuery = customQuery[:-4]

        # Examine subset of DataFrame of these duplicate energy runs
        tempdf = chan.query(customQuery)

        # Custom function to correctly compute the average yield (weighted by charge)
        avgY = lambda x : np.average(x, weights=tempdf.loc[x.index,'Yield err']**(-2))

        # Custom function to correctly compute the average of errors        #####CHECK THIS
        avgYErr = lambda x : np.average(np.ones_like(x), weights=tempdf.loc[x.index,'Yield err']**(-2))
        avgYErrEffcor = lambda x : np.average(np.ones_like(x), weights=tempdf.loc[x.index,'Yield err effcor']**(-2))

        # # TESTING STUFF ####
        # print(tempdf)
        # test1 = lambda x : np.average(x, weights=tempdf.loc[x.index,'Q int'])                       # WORKING
        # test2 = lambda x : np.sum( (tempdf.loc[x.index,'Q int']*tempdf.loc[x.index,'Yield err']) )  # WORKING
        # test3 = lambda x : np.sum( (tempdf.loc[x.index,'Q int']*tempdf.loc[x.index,'Yield err'])**2 )   # WORKING
        # test4 = lambda x : np.sqrt( np.sum( (tempdf.loc[x.index,'Q int']*     \
        #                         tempdf.loc[x.index,'Yield err'])**2 ) ) /     \
        #                         np.sum(tempdf.loc[x.index,'Q int'])   # WORKING
        #
        # # print(tempdf.groupby('Detector').agg({'Yield':test1}))
        # # print(tempdf.groupby('Detector').agg({'Yield':test2}))
        # # print(tempdf.groupby('Detector').agg({'Yield':test3}))
        # print(tempdf.groupby('Detector').agg({'Yield':test4}))
        # exit()


        # Custom function to give a new run number starting from 2000
        runnum = lambda x : newRunNum

        # Dictionary containing aggregation functions
        tempdict = {'Yield':avgY,'Yield err':avgYErr,
        'Area':'mean','Area err':avgYErr,'sig1':'mean','X2NDF':'mean',
        'IsValid':'mean','Status':'mean','Q int':'mean','Ep':'first',
        'Run':runnum,'Angle':'first','Yield effcor':avgY,
        'Yield err effcor':avgYErrEffcor
        }

        # New DataFrame that groups quantities by detector. Calculate new stuff
        # (yield, area, etc) and keep some old quantities,
        tempdf = tempdf.groupby('Detector').agg(tempdict)
        tempdf.reset_index(inplace=True)
        # print(tempdf)


        # # See what 'groupby' is actually doing
        # grouped_df = tempdf.groupby('Detector')
        # for key, item in grouped_df:
        #     print(grouped_df.get_group(key), '\n\n')
        # exit()


        # Combine the new averaged DataFrame with the old
        chan = pd.concat([chan,tempdf],sort=True, axis=0)

        # Drop the duplicate energy runs
        for runs in dupeRunList:
            customQuery = 'Run != %d'%runs
            chan.query(customQuery,inplace=True)

        newRunNum+=1


    if verbose > 0 :
        print('\n\nChecking if error in yield is greater than the yield POST OTHER STUFF:')
        print(chan.query('Yield < `Yield err`'))
        print('Number of points lost: ',len(chan.query('Yield < `Yield err`')))
    # Drop points where error in yield is greater than the yield
    chan.query('Yield > `Yield err`',inplace=True)


    # Pre-Cleaning complete, drop columns (no longer needed)
    chan.drop(columns = ['IsValid','Status','Area','Area err','sig1','X2NDF'], inplace=True)
    end_time = time.time()
    print('\n\nCleaning process required: %f seconds'%(end_time - start_time))

    # Final DataFrame
    # Assign column as index
    # chan = chan.set_index(['Run','Ep'])
    chan = chan.set_index(['Run'])
    print('\n\nFormat of final DataFrame:')
    print(chan.head())
    print('Size:',len(chan),'rows')


    # Save
    chan.to_csv('Yields/%s/%sYields.csv'%(ch.upper(),ch), sep=',')
    chan.to_excel('Yields/%s/%sYields.xlsx'%(ch.upper(),ch))

    # Make different sheets for each detector (quicker to plot and identify problems)
    with pd.ExcelWriter('Yields/%s/%sYieldsPerDetectors.xlsx'%(ch.upper(),ch)) as writer:
        for det in range(13):
            cut = (chan['Detector'] == det)
            tempdf = chan[cut]
            tempdf.to_excel(writer,sheet_name='%s'%det)

print('\n\n\nDONE!\n\n')

print('oooooooo O    O     O     O     O  O   O  Oooooooo    ')
print('   O     O    O    O O    O O   O  O  O   O           ')
print('   O     OooooO   OoooO   O  O  O  OOO    OooooooO    ')
print('   O     O    O  O     O  O   O O  O  O          O    ')
print('   O     O    O  O     O  O     O  O   O  oooooooO    ')
