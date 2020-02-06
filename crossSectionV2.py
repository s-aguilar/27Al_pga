import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12


# print('Attempting to create required directories: ')
# try:
#     os.mkdir('rMatrix')
#     os.mkdir('crossSection')
#     os.mkdir('crossSection/P1')
#     os.mkdir('crossSection/P2')
#     os.mkdir('crossSection/A1')
#     os.mkdir('crossPlots')
#     os.mkdir('crossPlots/P1')
#     os.mkdir('crossPlots/P2')
#     os.mkdir('crossPlots/A1')
# except OSError:
#     print ("Directories already exist!")
# else:
#     print('DONE!')


# Make `True` or `False` if you want plots
plots = False

channels = ['a1','p1','p2']

# Angle (deg) for each detector, negative is beam left, positive is beam right
#   Detector  :   Angle
#       00    :    120
#       01    :    105
#            ...
#       11    :   -105
#       12    :   -120


for ch in channels:

    print('\n\nWORKING ON %s CHANNEL:'%ch)

    # Read in DataFrame
    df = pd.read_excel('Yields/%s/%sYields.xlsx'%(ch.upper(),ch),index_col=[0])
    df.reset_index(inplace=True)

    # Keep necessary columns
    df = df[['Run','Angle','Ep','Yield effcor','Yield err effcor']]

    # Organize dataframe by energy
    df = df.sort_values(by=['Ep'])

    # Absolute value of `Angle` column
    df = df.assign(Angle=pd.Series(df['Angle'].abs(),index=df.index).values)

    # Rewrite Angles 120->60, 105->75
    df.loc[df['Angle']==120 , 'Angle'] = 60
    df.loc[df['Angle']==105 , 'Angle'] = 75

    # Convert energy to MeV
    df = df.assign(Ep=pd.Series(df['Ep']/1000,index=df.index).values)


    print('Starting averaging process...')
    start_time = time.time()

    # Custom aggregation functions to correctly calculate yields and their errors
    avgYEffcor = lambda x : np.average(x, weights=df.loc[x.index,'Yield err effcor']**(-2))
    avgYErrEffcor = lambda x : 1/np.sqrt(np.sum(df.loc[x.index,'Yield err effcor']**(-2)))

    tempdict = {'Ep':'first','Yield effcor':avgYEffcor,'Yield err effcor':avgYErrEffcor}

    # Newly averaged DataFrame
    averaged_df = df.groupby(['Run','Angle']).agg(tempdict)
    averaged_df.reset_index(inplace=True)

    end_time = time.time()
    print('Averaging Process required: %f seconds'%(end_time - start_time))


    # Stuff for computing Cross-Section
    solidAngle = 4*np.pi
    barn_conv = 1.e-24                               # Conversion of cm^2 to barn
    thickness = 6.e-6                                # 6 (ug/cm^2)
    density = 26.981                                # 27Al density (g/mol)
    numOfTarget = thickness/density*6.022e23        # thickness / (g/mol) * N_a (units of cm^2)
    numOfTarget *= barn_conv                        # Units of barn

    # Compute Cross-Section of averaged data
    Cross = averaged_df['Yield effcor'].values / numOfTarget / solidAngle
    Cross_err = averaged_df['Yield err effcor'].values / numOfTarget / solidAngle

    # Append Cross-sections back into pandas DataFrame
    averaged_df.drop(columns={'Run','Yield effcor','Yield err effcor'},inplace=True)
    averaged_df = averaged_df.assign(Cross=pd.Series(Cross,index=averaged_df.index).values)
    averaged_df = averaged_df.assign(Cross_err=pd.Series(Cross_err,index=averaged_df.index).values)
    averaged_df = averaged_df.rename(columns={'Cross_err': 'Cross err'})

    # Organize dataframe by energy and then by Angle number
    averaged_df = averaged_df.sort_values(by=['Ep','Angle'])

    averaged_df.set_index(['Ep'],inplace=True)
    averaged_df.reset_index(inplace=True)

    # Save finalized DataFrame
    # Make different sheets for each angle (quicker to plot and identify problems)
    with pd.ExcelWriter('%sCrossPerAngle.xlsx'%ch) as writer:
        for ang in np.array([0,15,30,45,60,75,90]):
            tempdf = averaged_df.query('Angle == %s'%ang)
            tempdf.to_excel(writer,sheet_name='%s'%ang)
            if plots:
                # Make the Cross-Section plot of averaged data
                plt.clf()

                temp_Eproton = tempdf['Ep'].values
                temp_Angle = tempdf['Angle'].values
                temp_Cross = tempdf['Cross'].values
                temp_Cross_err = tempdf['Cross err'].values

                plt.errorbar(temp_Eproton,temp_Cross,yerr=temp_Cross_err,fmt='k.',markersize='4') # was 2
                # plt.yscale('log')
                # plt.ylim(1e-6,1e-1)
                plt.xlim(2.05,3.25)
                plt.xlabel('$E_{p}$ (MeV)')
                plt.ylabel('Differential Cross-Section (barns/sr)')
                plt.title('%s %s$^{\circ}$'%(ch,ang))
                plt.savefig('plots/averaged/finalCrossSection/%s/%s_%s.png'%(ch.upper(),ch,ang),dpi=300)
                plt.yscale('log')
                plt.savefig('plots/averaged/finalCrossSection/%s/log_%s_%s.png'%(ch.upper(),ch,ang),dpi=300)


    # Saving DataFrame to `.dat` file
    np.savetxt('rMatrix/27Al_rMatrix_%s.dat'%ch, averaged_df.values,fmt='%5.4f   %d\t%f   %f')

print('\n\n\nDONE!\n\n')
print('oooooooo O    O     O     O     O  O   O  Oooooooo    ')
print('   O     O    O    O O    O O   O  O  O   O           ')
print('   O     OooooO   OoooO   O  O  O  OOO    OooooooO    ')
print('   O     O    O  O     O  O   O O  O  O          O    ')
print('   O     O    O  O     O  O     O  O   O  oooooooO    ')
