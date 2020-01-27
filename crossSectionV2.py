import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12


# Averages out run by run the left and right detectors, handles when there is
# only a right or left detector by just skipping that point.
# PRETTY SLOW
def average(df,runs):
    def averager(x,y,angie):
        try: # Both left and right
            left = temp.loc[temp['Detector'] == y,'Yield effcor'].values[0]
            left_err = temp.loc[temp['Detector'] == y,'Yield err effcor'].values[0]
            right = temp.loc[temp['Detector'] == x,'Yield effcor'].values[0]
            right_err = temp.loc[temp['Detector'] == x,'Yield err effcor'].values[0]

            yieldList.append( (left+right)/2. )
            yield_errList.append( (left_err**2 + right_err**2) ** .5 )
            energyList.append(temp['Ep'].values[0]/1000.)    # convert keV to MeV
            angleList.append(angie)
        except:
            try: # Only right
                right = temp.loc[temp['Detector'] == x,'Yield effcor'].values[0]
                right_err = temp.loc[temp['Detector'] == x,'Yield err effcor'].values[0]
                yieldList.append( right )
                yield_errList.append( right_err)
                energyList.append(temp['Ep'].values[0]/1000.)    # convert keV to MeV
                angleList.append(angie)
            except:
                try: # Only left
                    left = temp.loc[temp['Detector'] == y,'Yield effcor'].values[0]
                    left_err = temp.loc[temp['Detector'] == y,'Yield err effcor'].values[0]
                    yieldList.append( left )
                    yield_errList.append( left_err)
                    energyList.append(temp['Ep'].values[0]/1000.)    # convert keV to MeV
                    angleList.append(angie)
                except: return




    # Loop through each run, averaging the left and right detectors
    for _ in runs[:]:
        # tempp = df.set_index('Run')
        temp = df.query('Run == %s'%_)
        # print(temp.head(13))
        # exit()
        # print(temp.head(3))
        # print(temp.query('(Detector == 0) or (Detector == 12)'))
        averager(6,6,0)
        averager(5,7,15)
        averager(4,8,30)
        averager(3,9,45)
        averager(0,12,60)
        averager(1,11,75)
        averager(2,10,90)



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


detectors = ['det0','det1','det2','det3','det4','det5','det6','det7','det8',
             'det9','det10','det11','det12']


channels = ['a1','p1','p2']


# Angle (deg) for each detector, negative is beam left, positive is beam right
#                 00  01 02 03 04 05 06 07  08  09  10   11   12
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])



for ch in channels:
    # Read in DataFrame
    df = pd.read_excel('Yields/%s/%sYields.xlsx'%(ch.upper(),ch),index_col=[0])

    # Organize dataframe by energy and then by detector number
    df = df.sort_values(by=['Ep','Detector'])

    # cut is now DEPRECATED
    cut_df = df

    cut_Yield = cut_df['Yield'].values
    cut_Yield_err = cut_df['Yield err'].values

    cut_Yield_effcor = cut_df['Yield effcor'].values
    cut_Yield_err_effcor = cut_df['Yield err effcor'].values
    cut_Eproton = cut_df['Ep'].values/1000.    # convert MeV to keV
    cut_angle = cut_df['Angle'].values


    # Compute Cross-Section
    solidAngle = 4*np.pi
    barn_conv = 1e-24                               # Conversion of cm^2 to barn
    thickness = 11/(1e6)                             # 5ug/cm^2  maybe 5.7 check this ###################################
    numOfTarget = thickness*(1/26.981)*6.022e23     # thickness * (mol/27 g) * N_a (units of cm^2)
    numOfTarget *= barn_conv                        # Units of barn

    cut_Cross = cut_Yield_effcor / numOfTarget / solidAngle
    cut_Cross_err = cut_Yield_err_effcor / numOfTarget / solidAngle


    # """
    # NOT AVERAGED
    # Look at Before and After effect of cuts on the raw yield, and plot the post
    # cut, efficiency corrected cross section per detector
    for det in range(13):

        # Clear any other figure
        plt.clf()

        # BEFORE and AFTER cuts on raw yield per detector
        maskDet = ((df['Detector']==det))
        cut_maskDet = ((cut_df['Detector']==det))

        Yield = df['Yield'].values
        Yield_err = df['Yield err'].values
        Eproton = df['Ep'].values/1000.    # convert MeV to keV

        plt.errorbar(Eproton[maskDet],Yield[maskDet],yerr=Yield_err[maskDet],fmt='b.',markersize='2',label='Before')
        plt.errorbar(cut_Eproton[cut_maskDet],cut_Yield[cut_maskDet],yerr=cut_Yield_err[cut_maskDet],fmt='k.',markersize='2',label='After')

        plt.yscale('log')
        plt.ylim(1e-6,1)
        plt.xlim(2.05,3.3)
        # plt.xlim(2.35,3.3)
        plt.xlabel('$E_{p}$ (MeV)',fontsize=14)
        plt.ylabel('Yield (arb units)', fontsize=14)
        #plt.title('$^{27}$Al($\\mathrm{p,\\alpha_{1}\\gamma }$)$^{24}$Mg\t%s$^{\circ}$'%angle[det],fontsize=20)
        plt.title('$^{27}$Al($\\mathrm{p,X}$)\t%s$^{\circ}$'%angle[det],fontsize=20)
        plt.legend()
        plt.savefig('plots/notAveraged/yield/%s/%s_%s.png'%(ch.upper(),ch,detectors[det]),dpi=100)

        # Clean up
        plt.clf()

        # Plot the post cut, efficiency corrected cross section per detector
        cut_maskDet = ((cut_df['Detector']==det))

        # plt.plot(cut_Eproton[maskDet],cut_Cross[maskDet])
        plt.errorbar(cut_Eproton[cut_maskDet],cut_Cross[cut_maskDet],yerr=cut_Cross_err[cut_maskDet],fmt='k.',markersize='2')

        # print(len(cut_Eproton[cut_maskDet]),len(set(cut_Eproton[cut_maskDet])),len(cut_Eproton[cut_maskDet].duplicated()))
        plt.yscale('log')
        plt.ylim(1e-6,1)
        plt.xlim(2.05,3.3)
        plt.xlabel('$E_{p}$ (MeV)')
        plt.ylabel('Differential Cross-Section (barns/sr)')
        #plt.title('$^{27}$Al($\\mathrm{p,\\gamma \\alpha_{1}}$)$^{24}$Mg\t%s$^{\circ}$'%angle[det])
        plt.title('$^{27}$Al($\\mathrm{p,\\gamma X}$)\t%s$^{\circ}$'%angle[det])
        plt.savefig('plots/notAveraged/cross/%s/%s_%s.png'%(ch.upper(),ch,detectors[det]),dpi=300)

        # # SAVE YIELDS TO EXCEL TO MANUALLY PRUNE BAD RUNS LATER
        # df = pd.DataFrame(data=Yield[maskDet],index=Eproton[maskDet],columns=['Yield'])
        # df = df.assign(Yield_err=pd.Series(Yield_err[maskDet],index=df.index).values)
        # df.to_csv('yieldFiles/%s/%s_%s.csv'%(ch.upper(),ch,detectors[det]))
        # df.to_excel('yieldFiles/%s/%s_%s.xlsx'%(ch.upper(),ch,detectors[det]))

    with open('rMatrix/27Al_rMatrix_%s_allAngles.dat'%ch,'w') as f:
        for loop in range(len(cut_Cross)):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(cut_Eproton[loop],cut_angle[loop],cut_Cross[loop],cut_Cross_err[loop])
            f.write(printOut)
    # """


    # Array of unique run numbers in preserved order
    runArr = pd.unique(cut_df.index.values)

    # Average yields of L-R detectors
    energyList = []
    yieldList = []
    yield_errList = []
    angleList = []

    # start_time = time.time()
    # average(cut_df,runArr) #### VERY INEFFICIENT FUNCTION CALL
    # averaged_df = pd.DataFrame(data={'Ep':energyList, 'Yield effcor':yieldList, 'Yield err effcor':yield_errList, 'Angle':angleList})
    # end_time = time.time()
    # print('DONE!\t Process required: %f seconds'%(end_time - start_time))


    start_time = time.time()
    average(cut_df,runArr) #### VERY INEFFICIENT FUNCTION CALL
    averaged_df = pd.DataFrame(data={'Ep':energyList, 'Yield effcor':yieldList, 'Yield err effcor':yield_errList, 'Angle':angleList})
    end_time = time.time()
    print('Averaging Process required: %f seconds'%(end_time - start_time))


    # averaged_df.set_index(['Ep'],inplace=True)
    # print(len(df),len(cut_df),len(averaged_df))

    # Finalized cut and averaged data
    final_Yield_effcor = averaged_df['Yield effcor'].values
    final_Yield_err_effcor = averaged_df['Yield err effcor'].values
    final_Eproton = averaged_df['Ep'].values    # keV
    final_Angle = averaged_df['Angle'].values


    # Compute Cross-Section of finalized cut and averaged data
    final_Cross = final_Yield_effcor / numOfTarget / solidAngle
    final_Cross_err = final_Yield_err_effcor / numOfTarget / solidAngle


    # Stick cross sections back into pandas DataFrame
    averaged_df.drop(columns={'Yield effcor','Yield err effcor'},inplace=True)
    averaged_df = averaged_df.assign(Cross=pd.Series(final_Cross,index=averaged_df.index).values)
    averaged_df = averaged_df.assign(Cross_err=pd.Series(final_Cross_err,index=averaged_df.index).values)
    averaged_df = averaged_df.rename(columns={'Cross_err': 'Cross err'})


    # # Save finalized DataFrame
    # # Make different sheets for each angle (quicker to plot and identify problems)
    # with pd.ExcelWriter('%sCrossPerAngle.xlsx') as writer:
    #     for ang in np.array([0,15,30,45,60,75,90]):
    #         tempdf = averaged_df.query('Angle == %s'%ang)
    #         tempdf.set_index(['Ep'],inplace=True)
    #         tempdf.to_excel(writer,sheet_name='%s'%ang)


    # Create AZURE2 inputs
    # Create new output file / overwrite existing with empty file
    f = open('rMatrix/27Al_rMatrix_%s_cleaned.dat'%ch,'w')
    f.close()

    # Save finalized DataFrame
    with pd.ExcelWriter('%sCrossPerAngle.xlsx'%ch) as writer:
        for ang in np.array([0,15,30,45,60,75,90]):

            # Make the Cross-Section plot of finalized cut and averaged data
            plt.clf()

            # Select out the angle
            tempdf = averaged_df.query('Angle == %s'%ang)
            # print(len(tempdf))

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
            plt.show()

            # Append to file
            with open('rMatrix/27Al_rMatrix_%s_cleaned.dat'%ch,'a') as f:
                for loop in range(len(temp_Cross)):
                    printOut= '%f \t %d \t %.8f \t %.8f \n' %(temp_Eproton[loop],temp_Angle[loop],temp_Cross[loop],temp_Cross_err[loop])
                    f.write(printOut)


            # Make different sheets for each angle (quicker to plot and identify problems)
            tempdf.set_index(['Ep'],inplace=True)
            tempdf.to_excel(writer,sheet_name='%s'%ang)

        # exit()



print('DONE!')
