import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# This function checks to see if the errror is larger than the measurement,
# if it is it just artificially assigns it an error
def check(x,x_err):
    if x <= 1e-8:
        print('we got a problem')
    if x_err >= x:
        print('uhoh')
        print(x,'\t',x_err)
        return x,x*.9
    #if x_err >= x: return x_err
    else: return x,x_err


# Averages out run by run the left and right detectors, handles when there is
# only a right or left detector by just skipping that point.
# PRETTY SLOW
def average(df,runs):
    def averager(x,y,angie):
        try:
            left = temp.loc[temp['Detector'] == x,'Yield effcor'].values[0]
            left_err = temp.loc[temp['Detector'] == x,'Yield err effcor'].values[0]
        except: return
        try:
            right = temp.loc[temp['Detector'] == y,'Yield effcor'].values[0]
            right_err = temp.loc[temp['Detector'] == y,'Yield err effcor'].values[0]
        except: return

        yieldList.append( (left+right)/2. )
        yield_errList.append( (left_err**2 + right_err**2) ** .5 )
        energyList.append(temp['Ep'].values[0]/1000.)    # convert to MeV
        angleList.append(angie)

    for _ in runs[:]:
        temp = df.query('Run == %s'%_)
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

# Read in the data into dataframes
# dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')
# df1 = pd.read_csv('Yields/P1/p1Yields.csv')
# df2 = pd.read_csv('Yields/P2/p2Yields.csv')
df3 = pd.read_excel('Yields/A1/a1Yields.xlsx',index_col=[0])
# print(df3.head())


# effp1 = dfeff['p1'].values
# effp2 = dfeff['p2'].values
# effa1 = dfeff['a1'].values
# Angle = dfeff['Angle'].values

# Angle (deg) for each detector, negative is beam left, positive is beam right
#                 00  01 02 03 04 05 06 07  10  11  12   13   14
angle = np.array([120,105,90,45,30,15,0,-15,-30,-45,-90,-105,-120])

AnglesList=['0','15','30','45','90','105','120']


thickness = 5/(1e6)                             # 5ug/cm^2  maybe 5.7 check this ###################################
numOfTarget = thickness*(1/26.981)*6.022e23     # thickness * (mol/27 g) * N_a

q_e = 1.602e-19
scale = 1e-8    # 10^-8 C/pulse
q_corr = scale/(q_e)
barn_conv = 1/(1e-24)
solidAngle = 4*np.pi


# Example how to drop all data of a specific Run
# df3.query('Run != 97')


# Organize dataframe by energy and then by detector number
df3 = df3.sort_values(by=['Ep','Detector'])


# Cut low stats data
mask3Fit = ((df3['Area'] > 200)) #########
cut_df3 = df3.query('Area > 200')


# Array of unique run numbers in preserved order
runArr = pd.unique(cut_df3.index.values)


# Average yields of L-R detectors
energyList = []
yieldList = []
yield_errList = []
angleList = []
average(cut_df3,runArr) #### VERY INEFFICIENT FUNCTION CALL
averaged_df3 = pd.DataFrame(data={'Ep':energyList, 'Yield effcor':yieldList, 'Yield err effcor':yield_errList, 'Angle':angleList})
averaged_df3.set_index(['Ep'],inplace=True)


# a1Yield = averaged_df3['Yield'].values/ q_corr
# a1Yield_err = averaged_df3['Yield err'].values/ q_corr

exit()




a1Yield_effcor = a1Yield / effa1[0:len(a1Yield)]
a1Yield_err_effcor = a1Yield_err / effa1[0:len(a1Yield)]

a1Cross = a1Yield_effcor / numOfTarget * barn_conv / solidAngle
a1Cross_err = a1Yield_err_effcor / numOfTarget* barn_conv / solidAngle

a1Fit = averaged_df3['IsValid'].values
a1Eproton = averaged_df3['Ep'].values/1000    # Convert keV to MeV

# IsValid == 1 -> Good Fit
# IsValid == 0 -> Bad Fit
#
# Mask for which the fit was bad
# mask3Fit = ((df3['Area'] > 250) & (df3['IsValid'] == 1))


# # Sort by energy, keeping others consistent! DEPRECATED
# ind = a1Eproton.argsort()
# a1Eproton = a1Eproton[ind]
# a1Cross = a1Cross[ind]
# a1Cross_err = a1Cross_err[ind]
# a1Yield = a1Yield[ind]
# a1Yield_err = a1Yield_err[ind]


for det in range(len(detectors)):

    # Clear any other figure
    plt.clf()

    maskDet = ((df3['Detector']==detectors[det]))
    print(len(df3))
    print(maskDet)
    exit()

    plt.errorbar(a1Eproton[maskDet],a1Cross[maskDet],yerr=a1Cross_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    # plt.xlim(1.9,3.3)
    plt.xlim(2.35,3.3)
    plt.xlabel('$E_{p}$ (MeV)',fontsize=14)
    plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
    plt.title('$^{27}$Al($\\mathrm{p,\\alpha_{1}\\gamma }$)$^{24}$Mg\t%s$^{\circ}$'%angle[det],fontsize=20)
    plt.savefig('yieldPlots/A1/a1_%s.png'%detectors[det],dpi=100)

    # Clean up
    plt.clf()

    maskDet = ((df3['Detector']==detectors[det]) & mask3Fit)
    plt.plot(a1Eproton[maskDet],a1Cross[maskDet])
    plt.errorbar(a1Eproton[maskDet],a1Cross[maskDet],yerr=a1Cross_err[maskDet],fmt='b.',markersize='2')
    plt.errorbar(a1Eproton[maskDet],a1Yield[maskDet],yerr=a1Yield_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Cross-Section (barns/sr)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma \\alpha_{1}}$)$^{24}$Mg\t%s$^{\circ}$'%angle[det])
    plt.savefig('crossPlots/A1/a1_%s.png'%detectors[det],dpi=300)

    # SAVE YIELDS TO EXCEL TO MANUALLY PRUNE BAD RUNS LATER
    df = pd.DataFrame(data=a1Yield[maskDet],index=a1Eproton[maskDet],columns=['Yield'])
    df = df.assign(Yield_err=pd.Series(a1Yield_err[maskDet],index=df.index).values)
    df.to_csv('yieldFiles/A1/a1_%s.csv'%detectors[det])
    df.to_excel('yieldFiles/A1/a1_%s.xlsx'%detectors[det])
exit()
with open("rMatrix/27Al_rMatrix_a1_allAngles.dat","w") as f:
    for loop in range(len(a1Cross)):
        printOut= '%f \t %d \t %.8f \t %.8f \n' %(a1Eproton[loop],Angle[loop],a1Cross[loop],a1Cross_err[loop])
        f.write(printOut)



# """
f = open("rMatrix/27Al_rMatrix_a1.dat","w")
f.close()
for ang in AnglesList:
    _a1Eproton = []
    _Angle = []
    _a1Cross = []
    _a1Cross_err = []

    # Average out over same angle
    for x in range(867):    # total of 867 runs
        _a1Eproton.append(a1Eproton[int(13*x)])
        if ang == '0':
            # print(a1Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            cross = a1Cross[x*13+6]
            err = a1Cross_err[x*13+6]
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 ) # inflate errorbar
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            cross = (a1Cross[x*13+5]+a1Cross[x*13+7])/2
            err = (a1Cross_err[x*13+5]**2+a1Cross_err[x*13+7]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            cross = (a1Cross[x*13+4]+a1Cross[x*13+8])/2
            err = (a1Cross_err[x*13+4]**2+a1Cross_err[x*13+8]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            cross = (a1Cross[x*13+3]+a1Cross[x*13+9])/2
            err = (a1Cross_err[x*13+3]**2+a1Cross_err[x*13+9]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10])) )/2 )
            cross = (a1Cross[x*13+2]+a1Cross[x*13+10])/2
            err = (a1Cross_err[x*13+2]**2+a1Cross_err[x*13+10]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            cross = (a1Cross[x*13+1]+a1Cross[x*13+11])/2
            err = (a1Cross_err[x*13+1]**2+a1Cross_err[x*13+11]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            cross = (a1Cross[x*13]+a1Cross[x*13+12])/2
            err = (a1Cross_err[x*13]**2+a1Cross_err[x*13+12]**2)**.5
            cross,err = check(cross,err)
            _a1Cross.append(cross)
            _a1Cross_err.append( (err**2+(.05*_a1Cross[-1])**2)**.5 )


    _a1Eproton = np.array(_a1Eproton)
    _Angle = np.array(_Angle)
    _a1Cross = np.array(_a1Cross)
    _a1Cross_err = np.array(_a1Cross_err)

    # Sort by energy, keeping others consistent!
    ind = _a1Eproton.argsort()
    _a1Eproton = _a1Eproton[ind]
    _Angle = _Angle[ind]
    _a1Cross = _a1Cross[ind]
    _a1Cross_err = _a1Cross_err[ind]

    # Make the Cross-Section plot
    plt.clf()
    # plt.plot(_a1Eproton,_a1Cross)
    plt.errorbar(_a1Eproton,_a1Cross,yerr=_a1Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('a1 %s$^{\circ}$'%ang)
    plt.savefig('crossSection/A1/a1_%s.png'%ang,dpi=300)
    plt.clf()


    with open("rMatrix/27Al_rMatrix_a1.dat","a") as f:
        for loop in range(867):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_a1Eproton[loop],_Angle[loop],_a1Cross[loop],_a1Cross_err[loop])
            f.write(printOut)








"""
# Extract the columns of the DataFrame as numpy arrays
p1Run = df1['Run'].values
p1Det = df1['Detector'].values

p1Yield = df1['Yield'].values / q_corr
p1Yield_err = df1['Yield err'].values / q_corr

# p1Yield_effcor = p1Yield / effp1
# p1Yield_err_effcor = p1Yield_err / effp1

p1Yield_effcor = p1Yield / effp1[0:len(p1Yield)]
p1Yield_err_effcor = p1Yield_err / effp1[0:len(p1Yield)]

p1Cross = p1Yield_effcor / numOfTarget * barn_conv / solidAngle
p1Cross_err = p1Yield_err_effcor / numOfTarget* barn_conv / solidAngle

p1Fit = df1['IsValid'].values
p1Eproton = df1['Ep'].values/1000    # Convert keV to MeV


# IsValid == 1 -> Good Fit
# IsValid == 0 -> Bad Fit
#
# Mask for which the fit was bad
mask1Fit = ((df1['Area'] > 800) & (df1['IsValid'] == 1))

# Sort by energy, keeping others consistent!
ind = p1Eproton.argsort()
p1Eproton = p1Eproton[ind]
p1Cross = p1Cross[ind]
p1Cross_err = p1Cross_err[ind]
p1Yield = p1Yield[ind]
p1Yield_err = p1Yield_err[ind]

for det in range(len(detectors)):

    # Clear any other figure
    plt.clf()

    maskDet = ((df1['Detector']==detectors[det]))
    maskDet = maskDet[ind]

    plt.errorbar(p1Eproton[maskDet],p1Yield[maskDet],yerr=p1Yield_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Cross-Section (barns)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma p_{1}}$)$^{27}$Al\t%s$^{\circ}$'%angle[det])
    plt.savefig('yieldPlots/P1/p1_%s.png'%detectors[det],dpi=100)

    # Clean up
    plt.clf()

    maskDet = ((df1['Detector']==detectors[det]) & mask1Fit)
    maskDet = maskDet[ind]

    plt.plot(p1Eproton[maskDet],p1Cross[maskDet])
    plt.errorbar(p1Eproton[maskDet],p1Cross[maskDet],yerr=p1Cross_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Cross-Section (barns)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma p_{1}}$)$^{27}$Al\t%s$^{\circ}$'%angle[det])
    plt.savefig('crossPlots/P1/p1_%s.png'%detectors[det],dpi=100)


    # SAVE YIELDS TO EXCEL TO MANUALLY PRUNE BAD RUNS LATER
    df = pd.DataFrame(data=p1Yield[maskDet],index=p1Eproton[maskDet],columns=['Yield'])
    df = df.assign(Yield_err=pd.Series(p1Yield_err[maskDet],index=df.index).values)
    df.to_csv('yieldFiles/P1/p1_%s.csv'%detectors[det])
    df.to_excel('yieldFiles/P1/p1_%s.xlsx'%detectors[det])

with open("rMatrix/27Al_rMatrix_p1_allAngles.dat","w") as f:
    for loop in range(len(p1Cross)):
        printOut= '%f \t %d \t %.8f \t %.8f \n' %(p1Eproton[loop],Angle[loop],p1Cross[loop],p1Cross_err[loop])
        f.write(printOut)


# test2 = []
f = open("rMatrix/27Al_rMatrix_p1.dat","w")
f.close()
for ang in AnglesList:
    _p1Eproton = []
    _Angle = []
    _p1Cross = []
    _p1Cross_err = []

    # Average out over same angle
    for x in range(867):    # total of 867 runs
        _p1Eproton.append(p1Eproton[int(13*x)])
        if ang == '0':
            # print(p1Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            cross = p1Cross[x*13+6]
            err = p1Cross_err[x*13+6]
            cross,err = check(cross,err)
            _p1Cross.append(cross)
            # test2.append(p1Cross[x*13+6])
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5  )   # inflating the error bars for 5% systematic uncertainty
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            cross = (p1Cross[x*13+5]+p1Cross[x*13+7])/2
            # _p1Cross_err.append( (p1Cross_err[x*13+5]+p1Cross_err[x*13+7])/2 )
            err = (p1Cross_err[x*13+5]**2+p1Cross_err[x*13+7]**2)**.5
            cross,err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            cross = (p1Cross[x*13+4]+p1Cross[x*13+8])/2
            err = (p1Cross_err[x*13+4]**2+p1Cross_err[x*13+8]**2)**.5
            cross, err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            cross = (p1Cross[x*13+3]+p1Cross[x*13+9])/2
            err = (p1Cross_err[x*13+3]**2+p1Cross_err[x*13+9]**2)**.5
            cross, err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10]))/2 ) )
            cross = (p1Cross[x*13+2]+p1Cross[x*13+10])/2
            err = (p1Cross_err[x*13+2]**2+p1Cross_err[x*13+10]**2)**.5
            cross,err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            cross = (p1Cross[x*13+1]+p1Cross[x*13+11])/2
            err = (p1Cross_err[x*13+1]**2+p1Cross_err[x*13+11]**2)**.5
            cross,err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            cross = (p1Cross[x*13]+p1Cross[x*13+12])/2
            err = (p1Cross_err[x*13]**2+p1Cross_err[x*13+12]**2)**.5
            cross,err = check(cross,err)
            _p1Cross.append(cross)
            _p1Cross_err.append( (err**2+(.05*_p1Cross[-1])**2)**.5 )



    _p1Eproton = np.array(_p1Eproton)
    _Angle = np.array(_Angle)
    _p1Cross = np.array(_p1Cross)
    _p1Cross_err = np.array(_p1Cross_err)


    # Sort by energy, keeping others consistent!
    ind = _p1Eproton.argsort()
    _p1Eproton = _p1Eproton[ind]
    _Angle = _Angle[ind]
    _p1Cross = _p1Cross[ind]
    _p1Cross_err = _p1Cross_err[ind]


    # Make the Cross-Section plot
    plt.clf()
    plt.errorbar(_p1Eproton,_p1Cross,yerr=_p1Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma p_{3}}$)$^{27}$Al\t%s$^{\circ}$'%angle[det])
    plt.savefig('crossSection/P1/p1_%s.png'%ang,dpi=300)
    plt.clf()


    with open("rMatrix/27Al_rMatrix_p1.dat","a") as f:
        for loop in range(867):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_p1Eproton[loop],_Angle[loop],_p1Cross[loop],_p1Cross_err[loop])
            f.write(printOut)

# test1 = set(p1Cross[(df1['Detector']=='det_h0-6')][:16])
# print(test1)
# print('\n',test2[:16])
# test2 = set(test2[:16])
# if (test1==test2): print("its the same!")
# else: print("fuck its different!")

# """

"""
# Extract the columns of the DataFrame as numpy arrays
p2Run = df2['Run'].values
p2Det = df2['Detector'].values

p2Yield = df2['Yield'].values / q_corr
p2Yield_err = df2['Yield err'].values / q_corr

p2Yield_effcor = p2Yield / effp2[0:len(p2Yield)]
p2Yield_err_effcor = p2Yield_err / effp2[0:len(p2Yield)]

p2Cross = p2Yield_effcor / numOfTarget * barn_conv / solidAngle
p2Cross_err = p2Yield_err_effcor / numOfTarget* barn_conv / solidAngle

p2Fit = df2['IsValid'].values
p2Eproton = df2['Ep'].values/1000    # Convert keV to MeV

# IsValid == 1 -> Good Fit
# IsValid == 0 -> Bad Fit
#
# Mask for which the fit was bad
mask2Fit = ((df2['IsValid'] == 1) & (df2['Area'] > 800))

# Sort by energy, keeping others consistent!
ind = p2Eproton.argsort()
p2Eproton = p2Eproton[ind]
p2Cross = p2Cross[ind]
p2Cross_err = p2Cross_err[ind]
p2Yield = p2Yield[ind]
p2Yield_err = p2Yield_err[ind]
for det in range(len(detectors)):

    # Clear any other figure
    plt.clf()

    maskDet = ((df2['Detector']==detectors[det]))
    maskDet = maskDet[ind]

    plt.errorbar(p2Eproton[maskDet],p2Yield[maskDet],yerr=p2Yield_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Cross-Section (barns)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma p_{2}}$)$^{27}$Al\t%s$^{\circ}$'%angle[det])
    plt.savefig('yieldPlots/P2/p2_%s.png'%detectors[det],dpi=100)

    # Clean up
    plt.clf()

    maskDet = ((df2['Detector']==detectors[det]) & mask2Fit)
    maskDet = maskDet[ind]

    plt.clf()
    plt.plot(p2Eproton[maskDet],p2Cross[maskDet])
    plt.errorbar(p2Eproton[maskDet],p2Cross[maskDet],yerr=p2Cross_err[maskDet],fmt='b.',markersize='2')
    plt.errorbar(p2Eproton[maskDet],p2Yield[maskDet],yerr=p2Yield_err[maskDet],fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Cross-Section (barns/sr)')
    plt.title('$^{27}$Al($\\mathrm{p,\\gamma p_{2}}$)$^{27}$Al\t%s$^{\circ}$'%angle[det])
    plt.savefig('crossPlots/P2/p2_%s.png'%detectors[det],dpi=300)

    # SAVE YIELDS TO EXCEL TO MANUALLY PRUNE BAD RUNS LATER
    df = pd.DataFrame(data=p2Yield[maskDet],index=p2Eproton[maskDet],columns=['Yield'])
    df = df.assign(Yield_err=pd.Series(p2Yield_err[maskDet],index=df.index).values)
    df.to_csv('yieldFiles/P2/p2_%s.csv'%detectors[det])
    df.to_excel('yieldFiles/P2/p2_%s.xlsx'%detectors[det])


with open("rMatrix/27Al_rMatrix_p2_allAngles.dat","w") as f:
    for loop in range(len(p2Cross)):
        printOut= '%f \t %d \t %.8f \t %.8f \n' %(p2Eproton[loop],Angle[loop],p2Cross[loop],p2Cross_err[loop])
        f.write(printOut)



f = open("rMatrix/27Al_rMatrix_p2.dat","w")
f.close()
for ang in AnglesList:
    _p2Eproton = []
    _Angle = []
    _p2Cross = []
    _p2Cross_err = []

    # Average out over same angle
    for x in range(867):    # total of 867 runs
        _p2Eproton.append(p2Eproton[int(13*x)])
        if ang == '0':
            # print(p2Cross[x*13+6])
            _Angle.append(Angle[x*13+6])
            cross = p2Cross[x*13+6]
            err = p2Cross_err[x*13+6]
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '15':
            _Angle.append( (abs(Angle[x*13+5])+abs(Angle[x*13+7]))/2 )
            cross = (p2Cross[x*13+5]+p2Cross[x*13+7])/2
            err = (p2Cross_err[x*13+5]**2+p2Cross_err[x*13+7]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '30':
            _Angle.append( int( (abs(Angle[x*13+4])+abs(Angle[x*13+8]))/2 ) )
            cross = (p2Cross[x*13+4]+p2Cross[x*13+8])/2
            err = (p2Cross_err[x*13+4]**2+p2Cross_err[x*13+8]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '45':
            _Angle.append( int( (abs(Angle[x*13+3])+abs(Angle[x*13+9]))/2 ) )
            cross = (p2Cross[x*13+3]+p2Cross[x*13+9])/2
            err = (p2Cross_err[x*13+3]**2+p2Cross_err[x*13+9]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '90':
            _Angle.append( int( (abs(Angle[x*13+2])+abs(Angle[x*13+10]))/2 ) )
            cross = (p2Cross[x*13+2]+p2Cross[x*13+10])/2
            err = (p2Cross_err[x*13+2]**2+p2Cross_err[x*13+10]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '105':
            # _Angle.append( int( (abs(Angle[x*13+1])+abs(Angle[x*13+11]))/2 ) )
            _Angle.append( int( 75 ) )
            cross = (p2Cross[x*13+1]+p2Cross[x*13+11])/2
            err = (p2Cross_err[x*13+1]**2+p2Cross_err[x*13+11]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )
        elif ang == '120':
            # _Angle.append( int( (abs(Angle[x*13])+abs(Angle[x*13+12]))/2 ) )
            _Angle.append( int( 60 ) )
            cross = (p2Cross[x*13]+p2Cross[x*13+12])/2
            err = (p2Cross_err[x*13]**2+p2Cross_err[x*13+12]**2)**.5
            cross,err = check(cross,err)
            _p2Cross.append(cross)
            _p2Cross_err.append( (err**2+(.05*_p2Cross[-1])**2)**.5 )


    _p2Eproton = np.array(_p2Eproton)
    _Angle = np.array(_Angle)
    _p2Cross = np.array(_p2Cross)
    _p2Cross_err = np.array(_p2Cross_err)

    # Sort by energy, keeping others consistent!
    ind = _p2Eproton.argsort()
    _p2Eproton = _p2Eproton[ind]
    _Angle = _Angle[ind]
    _p2Cross = _p2Cross[ind]
    _p2Cross_err = _p2Cross_err[ind]


    # Make the Cross-Section plot
    plt.clf()
    # plt.plot(_p2Eproton,_p2Cross)
    plt.errorbar(_p2Eproton,_p2Cross,yerr=_p2Cross_err,fmt='b.',markersize='2')
    plt.yscale('log')
    plt.ylim(1e-6,1)
    plt.xlim(1.9,3.3)
    plt.xlabel('$E_{p}$ (MeV)')
    plt.ylabel('Differential Cross-Section (barns/sr)')
    plt.title('p2 %s$^{\circ}$'%ang)
    plt.savefig('crossSection/P2/p2_%s.png'%ang,dpi=300)
    plt.clf()

    with open("rMatrix/27Al_rMatrix_p2.dat","a") as f:
        for loop in range(867):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(_p2Eproton[loop],_Angle[loop],_p2Cross[loop],_p2Cross_err[loop])
            f.write(printOut)
# """





# """

print('DONE!')
