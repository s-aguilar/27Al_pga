import os
import time
import openpyxl # conda install -c anaconda openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from AtomicMassTable import GetElement
from chi2 import chi2_det,chi2_mat
from scipy.interpolate import LSQUnivariateSpline
from scipy import special
from scipy import interpolate
from scipy import constants
from scipy import integrate
import pdb

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# Use to switch directory paths to use (for the main or TESTING)
import os
currentDir = os.path.realpath('')
parentDir = os.path.realpath('..')



"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir


"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir

"""
Reads in the cross-section data produced by 'crossSectionV2.py' file and
performs a fit using legendre polynomial on the angular distributions.

The 'plot_it' variable enables/disables the plots demonstrating the different
order legendre polynomial fits 'a' analytic fit

"""


# GetElement returns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]
_barnsTocm2 = 1E-24         # cm2
_c = 2.99792458E10          # cm / s
_u = 931.494/_c**2          # MeV / c^2
_NA = 6.02214E23            # particles / mol


# Read in the data into dataframe
colNames = ['Energy','Angle','Cross-section','Error']
a1 = pd.read_table('rMatrix/27Al_rMatrix_a1_cleaned.dat',names=colNames)

dict_channels = {'a1':a1}


# Plotting: OFF = False, ON = True
plot_it = False


# Perform the analysis over all the channels
channels = ['a1']
for ch in channels:

    print(ch,' channel:\n-Working on Legendre fitting')
    angle = np.array([0,15,30,45,60,75,90])

    chan = dict_channels[ch]

    # is in Lab energy, convert to center-of-mass
    energyCM_chan = chan['Energy'].values*(_m27Al/(_m27Al+_m1H))   # Now its in E_cm
    chan = chan.assign(E_CM=pd.Series(energyCM_chan,index=chan.index).values)
    energy_chan = chan['Energy'].values
    angle_chan = chan['Angle'].values
    cross_chan = chan['Cross-section'].values
    crossErr_chan = chan['Error'].values


    # List of legendre orders to iterate over
    legendre_order = [[0],[0,1],[0,1,2],[0,1,2,3],[0,1,2,3,4],[0,1,2,3,4,5]]
    color_order = ['b','g','r','c','m','y']


    # Lists for analytic solutions
    ord_00,ord_01,ord_02,ord_03,ord_04, ord_05 = [],[],[],[],[],[]
    ord_err_00,ord_err_01,ord_err_02,ord_err_03,ord_err_04,ord_err_05 = [],[],[],[],[],[]
    ord_00_x2,ord_01_x2,ord_02_x2,ord_03_x2,ord_04_x2,ord_05_x2 = [],[],[],[],[],[]
    ord_00_x2ndf,ord_01_x2ndf,ord_02_x2ndf,ord_03_x2ndf,ord_04_x2ndf,ord_05_x2ndf = [],[],[],[],[],[]
    ord_00_pVal,ord_01_pVal,ord_02_pVal,ord_03_pVal,ord_04_pVal,ord_05_pVal = [],[],[],[],[],[]


    # Dicts for analytic solutions
    dict_ord_a = {'0':ord_00,'1':ord_01,'2':ord_02,'3':ord_03,'4':ord_04,'5':ord_05}
    dict_ord_err_a = {'0':ord_err_00,'1':ord_err_01,'2':ord_err_02,'3':ord_err_03,'4':ord_err_04,'5':ord_err_05}
    dict_ord_x2_a = {'0':ord_00_x2,'1':ord_01_x2,'2':ord_02_x2,'3':ord_03_x2,'4':ord_04_x2,'5':ord_05_x2}
    dict_ord_x2ndf_a = {'0':ord_00_x2ndf,'1':ord_01_x2ndf,'2':ord_02_x2ndf,'3':ord_03_x2ndf,'4':ord_04_x2ndf,'5':ord_05_x2ndf}
    dict_ord_pVal_a = {'0':ord_00_pVal,'1':ord_01_pVal,'2':ord_02_pVal,'3':ord_03_pVal,'4':ord_04_pVal,'5':ord_05_pVal}


    # Passing list to a set does not preserve the original order (ascending) so
    # manually sort it back into a new list
    # temp_nrg = set()
    # energyList = np.array([x for x in energy_chan if not (x in temp_nrg or temp_nrg.add(x))],dtype=np.float64)
    # print(energyList)
    energyList = pd.unique(energy_chan)
    # print(min(energyList),max(energyList))
    # exit()

    start_time = time.time()
    # perform Legendre fit over the angular distribution at each energy
    for nrg in energyList:

        # Get subset of DataFrame at each energy
        customQuery = 'Energy == %f'%nrg
        temp_chan = chan.query(customQuery)

        masked_nrg_chan = temp_chan['Energy'].values
        # masked_nrgCM_chan = temp_chan['E_CM'].values
        masked_ang_chan = np.radians(temp_chan['Angle'].values)
        masked_cross_chan = temp_chan['Cross-section'].values
        masked_err_chan = temp_chan['Error'].values

        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED


        # Up to Legendre Polynomial "BLAH" in legendre_order[list] index starts from 0, can go up to 5
        leg_ord = 4


        # Loop through the number of terms used in Legendre expansion:
        # -- a0, a0+a1, ... , a0+a1+a2+a3+a4, a0+a1+a2+a3+a4+a5
        for ind in range(0,leg_ord+1):

            results = chi2_mat(masked_cross_chan,masked_err_chan,masked_ang_chan,ind)

            dict_ord_a[str(legendre_order[ind][-1])].append(results[0])
            dict_ord_err_a[str(legendre_order[ind][-1])].append(results[1])
            dict_ord_x2_a[str(legendre_order[ind][-1])].append(results[2])
            dict_ord_x2ndf_a[str(legendre_order[ind][-1])].append(results[3])
            dict_ord_pVal_a[str(legendre_order[ind][-1])].append(results[4])

            if plot_it:
                # Note: The legendre.legval function requires the coefficients
                # of all terms up to whatever order you want, must pad the even
                # terms with 0's, this is done through the 'convert' function
                temp_coef = results[0]
                temp_coef_err = results[1]

                # Recording points of some angular distributions
                # # if nrg == 4.30930:
                # if nrg < 4.962 and nrg > 4.96:
                #     print(ch,legendre_order[ind],nrg)
                #     # df = pd.DataFrame(data=np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),index=ang_knots,columns=['legVal'])
                #     df = pd.DataFrame(data=ang_knots,columns=['angle'])
                #     df = df.assign(legVal=pd.Series(np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),index=df.index).values)
                #     df.to_excel('legendre_out/DATA/analytically/%s/a%d/%s_a%d_legendrePointsAtSingleEnergy.xlsx'%(ch,legendre_order[ind][-1],ch,legendre_order[ind][-1]))

                # Scatter has bug when the y-values are small, the y-axis does not autoscale
                # plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])
                # if ind ==2:
                #     plt.plot(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],label='%d$^{th}$ Order'%legendre_order[ind][-1])
                plt.plot(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],label='Order - %d'%legendre_order[ind][-1])
        # Angular distribution plots
        if plot_it:
            plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
            plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
            plt.xlabel('Angle',fontsize=14)
            plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
            plt.title('$^{24}$Mg($\\alpha$,p$_{%s}\\gamma$) - %5.4f MeV' % (ch[-1],nrg),fontsize=20)
            plt.legend(prop={'size': 14})
            plt.savefig('legendre_out/fits/%s/%5.4fMeVFit.png'%(ch,nrg),dpi=100)
            plt.clf()
    end_time = time.time()
    print('Process required: %f seconds'%(end_time - start_time))

    # Convert lists into arrays
    for ind in range(0,leg_ord+1):
        dict_ord_a[str(legendre_order[ind][-1])] = np.array(dict_ord_a[str(legendre_order[ind][-1])])
        dict_ord_err_a[str(legendre_order[ind][-1])] = np.array(dict_ord_err_a[str(legendre_order[ind][-1])])


    # Now we have all the coefficients for each order fit, lets plot them
    # plot 'a' coefficients for each order as function of energy as well as
    # the chi2ndf
    ord_step = 0
    print('-Working on plotting the coefficients')
    start_time = time.time()
    for ind in range(0,leg_ord+1):
        for ord_step in range(0,leg_ord+1):
            plt.clf()
            if ind > ord_step: continue
            a = [ (_[ind],) for _ in dict_ord_a[str(legendre_order[ord_step][-1])]]     # analytic solution
            plt.scatter(energyList,a,s=8)
            plt.title('%s channel - Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1],legendre_order[ord_step][-1]))
            plt.savefig('legendre_out/coef_curve/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=100)
            ord_step += 1

        plt.clf()

        x2 = dict_ord_x2_a[str(legendre_order[ind][-1])]
        x2ndf = dict_ord_x2ndf_a[str(legendre_order[ind][-1])]
        # print(str(legendre_order[ind][-1]))
        plt.scatter(energyList,x2ndf,s=8)
        plt.title('%s channel - $\chi^{2}$ for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1]))
        plt.savefig('legendre_out/chi2/%s/a%d/chi2.png'%(ch,legendre_order[ind][-1]),dpi=100)

    end_time = time.time()
    print('Process required: %f seconds'%(end_time - start_time))

    plt.clf()

    # """
    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    start_time = time.time()
    colLabels = ['a0','a1','a2','a3','a4','a5']
    colLabelss = ['a0_err','a1_err','a2_err','a3_err','a4_err','a5_err']
    for _ in range(0,leg_ord+1):
        stuff = dict_ord_a[str(legendre_order[_][-1])]
        stuff2 = dict_ord_err_a[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_+1])
        df1 =  pd.DataFrame(data=stuff2,index=energyList,columns=colLabelss[0:_+1])
        df = pd.concat([df, df1], axis=1, sort=False)
        df = df.assign(Chi2=pd.Series(dict_ord_x2_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Pval=pd.Series(dict_ord_pVal_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Chi2NDF=pd.Series(dict_ord_x2ndf_a[str(legendre_order[_][-1])],index=df.index).values)
        df.to_csv('legendre_out/DATA/%s/a%d/%s_a%dFit.csv'%(ch,legendre_order[_][-1],ch,legendre_order[_][-1]))
        df.to_excel('legendre_out/DATA/%s/a%d/%s_a%dFit.xlsx'%(ch,legendre_order[_][-1],ch,legendre_order[_][-1]))

    end_time = time.time()
    print('Process required: %f seconds'%(end_time - start_time))
    # """

    start_time = time.time()
    # Now I have all the 'a_i' coefficients for all the many order fits, need
    # to select/filter the necessary 'a0' terms from the fits. This will be
    # done through the p-value.
    # Read in all my a0s and p-values and compare each energies fit side by side
    # starting from lowest order to highest, select order which pvalue is maximal,
    # append the a0 term + continue and repeat for the next energy step
    # (NOT USED ANYMORE)

    a0_final = []
    a0_err_final = []
    p_final = []
    x2_final = []
    order = []

    # Only keep the a0s from the 1st order legendre polynomial fit
    for ind in range(len(ord_00)):
        a0_final.append(dict_ord_a[str(legendre_order[1][-1])][ind][0])
        a0_err_final.append(dict_ord_err_a[str(legendre_order[1][-1])][ind][0])
        p_final.append(dict_ord_pVal_a[str(legendre_order[1][-1])][ind])
        order.append(str(legendre_order[1][-1]))


    # Save the points
    df = pd.DataFrame(data=a0_final,index=energyList)
    df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).values)
    df = df.assign(Pval=pd.Series(p_final,index=df.index).values)
    df = df.assign(Order=pd.Series(order,index=df.index).values)
    df.to_excel('legendre_out/DATA/%s/%sFits.xlsx'%(ch,ch))

    end_time = time.time()
    print('Process required: %f seconds'%(end_time - start_time))




    # Save angle integrated cross section for rMatrix
    with open('rMatrix/27Al_rMatrix_a1_angInt.dat','w') as f:
        print(len(a0_final),len(a0_err_final),len(energyList))

        cross = np.array(a0_final)*4*np.pi
        cross_err = np.array(a0_err_final)*4*np.pi
        print(len(a0_final),len(a0_err_final),len(energyList))
        # exit()
        for loop in range(len(energyList)):
            printOut= '%f \t %d \t %.8f \t %.8f \n' %(energyList[loop],0,cross[loop],cross_err[loop])
            f.write(printOut)





    exit()






















    # Now spline the 'a0' coefficients and plot as function of energy and overlay
    # the original 'a0' coefficents and make sure the splining is done well.

    a0_final = np.array(a0_final)
    a0_err_final = np.array(a0_err_final)

    solidAngle = 4*np.pi

    # THIS IS FOR LAB ENERGIES
    # I want to spline from # <= 4.03514 and spline # > 4.26996 THIS IS FOR LAB ENERGY
    # e_knots_mask1 = energyList <= 4.03514
    # spline1_mask = energyList <= 4.03514
    # e_knots_mask2 = energyList >= 4.26996
    # spline2_mask = energyList >= 4.26996

    # FOR CM ENERGIES
    # e_knots_mask1 = energyList <= 3.459
    # spline1_mask = energyList <= 3.459
    # e_knots_mask2 = energyList >= 3.659
    # spline2_mask = energyList >= 3.659

    # e_mask1 = energyList <= 3.459
    # e_mask2 = energyList >= 3.65


    if ch == 'p1':
        # e_mask1 = energyList >= 3.745
        e_mask1 = energyList >= 3.745


        interior1 = energyList[e_mask1]
        e1spline_x = np.linspace(min(interior1),max(interior1),1000)

        # Must optimize interior knots 't' such as not to produce spurious curves
        # knots1 = interior1[3:-1:2]
        # knots1 = np.ma.concatenate([interior1[3:24:2],interior1[25:-1:2]])
        # knots1 = np.ma.concatenate([interior1[3:24:2],interior1[25:30:1],interior1[31:-1:2]])

        knots1 = np.ma.concatenate([interior1[3:24:2],interior1[25:30:1],interior1[32:37:2],interior1[37:-1:2]])


        weights1 = 1/(a0_err_final[e_mask1]*solidAngle)            # 1/sigma

        # Create a cubic spline representation that is weighted by the uncertainties 'w'
        spline1 = LSQUnivariateSpline(energyList[e_mask1],a0_final[e_mask1]*solidAngle, t = knots1, w = weights1)
        y_val1 = spline1(e1spline_x)

    elif ch == 'p2':
        e_mask1 = energyList <= 3.459
        e_mask2 = energyList >= 3.67
        interior1 = energyList[e_mask1]
        e1spline_x = np.linspace(min(interior1),max(interior1),100)

        # Must optimize interior knots 't' such as not to produce spurious curves
        knots1 = interior1[2:-1:2]
        weights1 = 1/(a0_err_final[e_mask1]*solidAngle)            # 1/sigma

        # Create a cubic spline representation that is weighted by the uncertainties 'w'
        spline1 = LSQUnivariateSpline(energyList[e_mask1],a0_final[e_mask1]*solidAngle, t = knots1, w = weights1)
        y_val1 = spline1(e1spline_x)


        interior2 = energyList[e_mask2]
        e2spline_x = np.linspace(min(interior2),max(interior2),1000)

        # Must optimize interior knots 't' such as not to produce spurious curves
        # knots2 = interior2[1:-1:2]
        knots2 = np.ma.concatenate([interior2[1:6:3],interior2[7:18:1],interior2[19:46:2],interior2[46:65:1],interior2[67:84:2],interior2[85:93:3],interior2[94:-1:2]])
        # knots2 = np.ma.concatenate([interior2[1:6:3],interior2[7:18:1],interior2[19:46:2],interior2[46:51:2],interior2[51:54:2],interior2[55:65:1],interior2[67:84:2],interior2[85:93:3],interior2[94:-1:2]])
        weights2 = 1/(a0_err_final[e_mask2]*solidAngle)            # 1/sigma

        # Create a cubic spline representation that is weighted by the uncertainties 'w'
        spline2 = LSQUnivariateSpline(energyList[e_mask2],a0_final[e_mask2]*solidAngle, t = knots2, w = weights2)
        y_val2 = spline2(e2spline_x)

    # Plot the splines and the points
    plt.clf()

    plt.plot(e1spline_x,y_val1,c='b',label='Cubic Spline')
    plt.errorbar(energyList[e_mask1],a0_final[e_mask1]*solidAngle,yerr=a0_err_final[e_mask1]*solidAngle,fmt='.',c='k')
    if ch == 'p2':
        plt.plot(e2spline_x,y_val2,c='b')
        plt.errorbar(energyList[e_mask2],a0_final[e_mask2]*solidAngle,yerr=a0_err_final[e_mask2]*solidAngle,fmt='.',c='k')

    plt.xlabel('E$_{CM}$ (MeV)',fontsize=14)
    plt.ylabel('Angle Integrated Cross-Section (barns)', fontsize=14)
    plt.yscale('log')
    plt.legend()
    plt.title('%s channel - a0*4$\pi$ vs E - Excitation Curve' % ch)
    plt.savefig('legendre_out/excitationCurve/%s/%s_a0Curve.png' % (ch,ch),dpi=300)
    # plt.show()
    plt.clf()
    # exit()


    # Save the points
    df = pd.DataFrame(data=a0_final,index=energyList,columns=['a0'])
    df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).values)
    df = df.assign(fitOrder=pd.Series(order,index=df.index).values)
    df.to_excel('legendre_out/DATA/%s/%sCrossPoints.xlsx'%(ch,ch))

    # Splined CrossSection, first merge the splines
    if ch == 'p1':
        yvals = y_val1
        xvals = e1spline_x
    else:
        yvals = np.ma.concatenate([y_val1,y_val2])
        xvals = np.ma.concatenate([e1spline_x,e2spline_x])
    df = pd.DataFrame(data=xvals,columns=['splineX'])
    df = df.assign(splineY=pd.Series(yvals,index=df.index).values)
    df.to_excel('legendre_out/DATA/%s/%sSplinePoints.xlsx'%(ch,ch))



    # Calculate the reaction rates and write them as both 'csv' and 'xlsx' format
    print('-Calculating and writing reaction rates into .csv and .xlsx files')
    fname = 'iliadis_Temps.dat'
    # fname = 'rate_Temps.dat'
    temperature = np.loadtxt( fname )

    # reduced mass
    mu = _m24Mg*_m4He/(_m24Mg+_m4He)

    rxnRate = []
    integrandList = []
    eList = []
    filenames = []
    tempNames= []

    # temp variable contains the list of integrands and energies which will be written to an excel file
    for T in temperature:

        plt.clf()

        E = energyList
        integrand = (a0_final*solidAngle) * E * np.exp(-11.604*E/T)
        temp = integrand[e_mask1]
        integrand_err = (a0_err_final*solidAngle) * E * np.exp(-11.604*E/T) # MAYBE?


        # Combine integral of both splines
        E1 = energyList[e_mask1]
        tempE = E1
        weights1 = 1/integrand_err[e_mask1]            # 1/sigma

        # Create a cubic spline representation that is weighted by the uncertainties 'w'
        spline1 = LSQUnivariateSpline(E1,integrand[e_mask1], t = knots1, w = weights1)
        y_val1 = spline1(e1spline_x)

        if ch == 'p2':
            E2 = energyList[e_mask2]
            weights2 = 1/integrand_err[e_mask2]            # 1/sigma

            temp = np.ma.concatenate([integrand[e_mask1],integrand[e_mask2]])
            tempE = np.ma.concatenate([E1,E2])

            # Create a cubic spline representation that is weighted by the uncertainties 'w'
            spline2 = LSQUnivariateSpline(E2,integrand[e_mask2], t = knots2, w = weights2)
            y_val2 = spline2(e2spline_x)
            plt.plot(e2spline_x,y_val2,c='b')


        # Stuff for integrand
        tempNames.append('T=%3.2fK'%T)
        integrandList.append(temp)
        eList.append(tempE)



        plt.plot(e1spline_x,y_val1,c='b',label='Cubic Spline')

        # plt.xlim(0,10)
        plt.errorbar(E[e_mask1],integrand[e_mask1],yerr=integrand_err[e_mask1],fmt='.',c='k')
        if ch == 'p2':
            plt.errorbar(E[e_mask2],integrand[e_mask2],yerr=integrand_err[e_mask2],fmt='.',c='k')
            plt.plot(e2spline_x,y_val2,c='b',label='Cubic Spline')
        plt.yscale('log')
        plt.ylabel('Reaction Rate Integrand',fontsize=14)
        plt.xlabel('E$_{CM}$ (MeV)',fontsize=14)
        plt.title('%s - Integrand at %fGK'%(ch,T),fontsize=20)
        plt.legend()
        filenames.append('%s_DataReactionRate_integrand_at_%fGK.png'%(ch,T))
        plt.savefig('%s_DataReactionRate_integrand_at_%fGK.png'%(ch,T),dpi=300)
        # plt.show()

        if ch == 'p1':
            yvals = y_val1
            xvals = e1spline_x
        else:
            yvals = np.ma.concatenate([y_val1,y_val2])
            xvals = np.ma.concatenate([e1spline_x,e2spline_x])
        df = pd.DataFrame(data=xvals,columns=['splineX'])
        df = df.assign(splineY=pd.Series(yvals,index=df.index).values)
        df.to_excel('legendre_out/DATA/%s/%sIntegrandSplinePoints%f.xlsx'%(ch,ch,T))


        integrand_x = E[e_mask1]
        integrand_y = integrand[e_mask1]
        if ch == 'p2':
            integrand_x = np.ma.concatenate( [integrand_x,E[e_mask2]] )
            integrand_y = np.ma.concatenate( [integrand_y,integrand[e_mask2]] )

        df = pd.DataFrame(data=integrand_x,columns=['integrand_x'])
        df = df.assign(integrand_y=pd.Series(integrand_y,index=df.index).values)
        df.to_excel('legendre_out/DATA/%s/%sIntegrandDataPoints%f.xlsx'%(ch,ch,T))


        intg = spline1.integral(min(E1),max(E1))
        if ch == 'p2':
            intg += spline2.integral(min(E2),max(E2))


        # rate =  3.7313E10 * mu**(-.5) * T**(-1.5) * intg
        rate =  rxnRateCONST * mu**(-.5) * T**(-1.5) * intg
        rxnRate.append(rate)

    # Make different sheets for integrand for each temperature
    with pd.ExcelWriter('legendre_out/integrand/%s_integrands.xlsx'%ch) as writer:
        for ind in range(len(tempNames)):
            print(tempNames[ind])
            tempdf = pd.DataFrame(data=integrandList[ind],index=eList[ind],columns=['Integrand'])
            tempdf = tempdf.rename_axis('E_cm (MeV)',axis=0) # rename the index axis
            print(tempdf.head(5))
            tempdf.to_excel(writer,sheet_name=tempNames[ind])


    import imageio
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave('legendre_out/%s_data_integrand.gif'%ch, images,duration=1)
    # exit()

    plt.clf()
    rxnRate = np.array(rxnRate)
    plt.plot(temperature,rxnRate)
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-20,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('%s - Reaction Rate'%ch,fontsize=20)
    plt.savefig('%s_DataReactionRate_integrate.png'%ch,dpi=300)
    plt.clf()
    df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)
    df.to_csv('legendre_out/DATA/%s/a%d/%s_rates.csv'%(ch,legendre_order[0][-1],ch))
    df.to_excel('legendre_out/DATA/%s/a%d/%s_rates.xlsx'%(ch,legendre_order[0][-1],ch))

    # exit()

    # continue
    print('\n\n')
# """

print('DONE!')
