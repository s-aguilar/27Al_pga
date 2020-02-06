'''
Script reads in my 27Al data and Nelson data, and plots them. Noticed an approx.
uniform energy shift between our data so I shift mine +4 keV.

Resonances thats widths were guessed by eye are now properly fit to extract a
more quantitative value for the widths. These widths are then used to calculate
a more accurate and precise 27Al target we used for our data
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussFit(x,x0,A,sig):
    """
    Normalized gaussian that is used to fit resonance:
        x   - position
        x0  - centroid position
        A   - peak area
        sig - width
    """
    coef = A / (sig*np.sqrt(2*np.pi))
    gaussPart = np.exp(-.5*((x-x0)/sig)**2)

    return coef * gaussPart

def gaussFitPlusBack(x,x0,A,sig,c0,c1):
    """
    Normalized gaussian that is used to fit resonance:
        x   - position
        x0  - centroid position
        A   - peak area
        sig - width

    Linear background:
        c0 + c1*x
    """
    coef = A / (sig*np.sqrt(2*np.pi))
    gaussPart = np.exp(-.5*((x-x0)/sig)**2)
    back = c0 + c1*x

    return (coef * gaussPart) + back

def lorenFit(x,x0,gam,A):
    """
    Pseudo-Lorentz/Cauchy distribution
        x   - position
        x0  - centroid position
        gam - HWHM (half width half max)
        A   - peak area
    """
    coef = A*.5*gam/np.pi
    lorenPart = ( (x-x0)**2 + .25*gam**2)**(-1)

    return coef * lorenPart

def lorenFitPlusBack(x,x0,gam,A,c0,c1):
    """
    Pseudo-Lorentz/Cauchy distribution
        x   - position
        x0  - centroid position
        gam - HWHM (half width half max)
        A   - peak area

    Linear background:
        c0 + c1*x
    """
    coef = A*.5*gam/np.pi
    lorenPart = ( (x-x0)**2 + .25*gam**2)**(-1)
    back = c0 + c1*x

    return (coef * lorenPart) + back



# Read in mine and Nelson data into DataFrames
colNames = ['Energy','Angle','Cross-section','Error']
me_a1 = pd.read_table('rMatrix/27Al_rMatrix_a1_angInt.dat',names=colNames)
nel_a1 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)

# Nelson a1 angles: 120, 135, 150
nel_a1_120 = nel_a1.query('Angle == 120')
nel_a1_135 = nel_a1.query('Angle == 135')
nel_a1_150 = nel_a1.query('Angle == 150')

# My data with a +4 keV shift
me_a1_4shift = me_a1.apply(lambda x: x+.004 if x.name == 'Energy' else x)


# plt.errorbar(me_a1['Energy'],me_a1['Cross-section'],yerr=me_a1['Error'],fmt='.',c='r',label='Me')
# plt.errorbar(me_a1_4shift['Energy'],me_a1_4shift['Cross-section'],yerr=me_a1_4shift['Error'],fmt='.',c='b',label='Me shift')
# plt.errorbar(nel_a1_120['Energy'],nel_a1_120['Cross-section'],yerr=nel_a1_120['Error'],fmt='.',c='k',label='Nel120')
# plt.errorbar(nel_a1_135['Energy'],nel_a1_135['Cross-section'],yerr=nel_a1_135['Error'],fmt='.',c='g',label='Nel135')
# plt.errorbar(nel_a1_150['Energy'],nel_a1_150['Cross-section'],yerr=nel_a1_150['Error'],fmt='.',c='c',label='Nel150')
# plt.xlim(2.295,2.95)
# plt.yscale('log')
# plt.legend()
# plt.show()
# exit()

# """
# Weighted fit of my a1 data
me_a1_4shift_slice = me_a1_4shift.query('Energy >= 2.735 and Energy <= 2.745')
# print(me_a1_4shift_slice['Error']/me_a1_4shift_slice['Cross-section'])
print(me_a1_4shift_slice['Error'].values/me_a1_4shift_slice['Cross-section'].values)
# print(me_a1_4shift_slice['Error'].values)
# print(me_a1_4shift_slice['Cross-section'].values)
# exit()

# Gauss Fit
ePoints = np.linspace(min(me_a1_4shift_slice['Energy']),max(me_a1_4shift_slice['Energy']),1000)
p0 = np.array([2.7405,1e-4,0.001]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=gaussFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = gaussFit(ePoints, *popt1)    # Fit line

p00 = np.array([2.7405,1e-4,0.001,1,1]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=gaussFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = gaussFitPlusBack(ePoints, *popt11)    # Fit line

# print('My G fit:',popt1)
# print('My G+B fit:',popt11)
plt.errorbar(me_a1_4shift_slice['Energy'],me_a1_4shift_slice['Cross-section'],
            yerr=me_a1_4shift_slice['Error'],fmt='.',c='k',label='Data')
plt.plot(ePoints,yfit1,c='r',label='G')
plt.plot(ePoints,yfit11,c='b',label='G+B')
# plt.yscale('log')
plt.text(2.735,0.043,'Pub. Width: 0.6 keV')
# print(type(popt1[1]))
# exit()
w1 = 2.355*popt1[2]*1000
w11 = 2.355*popt11[2]*1000
# print('My G width:',w1)
# print('My G+B width:',w11)
plt.text(2.735,0.04,'G Width: %4.3f keV'%w1)
plt.text(2.735,0.037,'G+B Width: %4.3f keV'%w11)
plt.title('Fit of an a$_{1}$ Resonance in $^{27}$Al Data')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Cross-section (arb. units)')
# plt.legend()
# plt.show()


# Lorentz Fit
p0 = np.array([2.7405,1e-4,1]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=lorenFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = lorenFit(ePoints, *popt1)    # Fit line

p00 = np.array([2.7405,1e-4,1,1,1]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=lorenFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = lorenFitPlusBack(ePoints, *popt11)    # Fit line

# print('My L fit:',popt1)
# print('My L+B fit:',popt11)
plt.plot(ePoints,yfit1,c='g',label='L')
plt.plot(ePoints,yfit11,c='m',label='L+B')
w1 = popt1[1]*1000
w11 = popt11[1]*1000
print('My L width:',w1)
print('My L+B width:',w11)
plt.text(2.735,0.034,'L Width: %4.3f keV'%w1)
plt.text(2.735,0.031,'L+B Width: %4.3f keV'%w11)
plt.legend()
plt.xlim(2.7345,2.745)
# plt.yscale('log')
plt.savefig('resWidthAnal.png',dpi=900)
# plt.show()
plt.clf()



# Weighted fit of Nelson a1 data
me_a1_4shift_slice = nel_a1_120.query('Energy >= 2.735 and Energy <= 2.745')

# Gauss Fit
ePoints = np.linspace(min(me_a1_4shift_slice['Energy']),max(me_a1_4shift_slice['Energy']),1000)
p0 = np.array([2.7405,1e-4,0.001]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=gaussFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = gaussFit(ePoints, *popt1)    # Fit line

p00 = np.array([2.7405,1e-4,0.001,1,1]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=gaussFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = gaussFitPlusBack(ePoints, *popt11)    # Fit line

# print('My G fit:',popt1)
# print('My G+B fit:',popt11)
plt.errorbar(me_a1_4shift_slice['Energy'],me_a1_4shift_slice['Cross-section'],
            yerr=me_a1_4shift_slice['Error'],fmt='.',c='k',label='Data')
plt.plot(ePoints,yfit1,c='r',label='G')
plt.plot(ePoints,yfit11,c='b',label='G+B')
# plt.yscale('log')
plt.text(2.735,0.007,'Pub. Width: 0.6 keV')
# print(type(popt1[1]))
# exit()
w1 = 2.355*popt1[2]*1000
w11 = 2.355*popt11[2]*1000
# print('My G width:',w1)
# print('My G+B width:',w11)
plt.text(2.735,0.0065,'G Width: %4.3f keV'%w1)
plt.text(2.735,0.006,'G+B Width: %4.3f keV'%w11)
plt.title('Fit of an a$_{1}$ Resonance in $^{27}$Al Nelson Data 120$^{o}$')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Cross-section (arb. units)')
# plt.legend()
# plt.show()


# Lorentz Fit
p0 = np.array([2.7405,1e-4,1]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=lorenFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = lorenFit(ePoints, *popt1)    # Fit line

p00 = np.array([2.7405,1e-4,1,1,1]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=lorenFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = lorenFitPlusBack(ePoints, *popt11)    # Fit line

# print('My L fit:',popt1)
# print('My L+B fit:',popt11)
plt.plot(ePoints,yfit1,c='g',label='L')
plt.plot(ePoints,yfit11,c='m',label='L+B')
w1 = popt1[1]*1000
w11 = popt11[1]*1000
# print('My L width:',w1)
# print('My L+B width:',w11)
plt.text(2.735,0.0055,'L Width: %4.3f keV'%w1)
plt.text(2.735,0.005,'L+B Width: %4.3f keV'%w11)
plt.xlim(2.7345,2.745)
plt.legend()
# plt.yscale('log')
plt.savefig('resWidthAnalNelson120.png',dpi=900)
# plt.show()
plt.clf()
# """


























# Other resonance
print('\n\nNEW RES')

# Weighted fit of my a1 data
me_a1_4shift_slice = me_a1_4shift.query('Energy >= 2.705 and Energy <= 2.72')


# Gauss Fit
ePoints = np.linspace(min(me_a1_4shift_slice['Energy']),max(me_a1_4shift_slice['Energy']),1000)
p0 = np.array([2.71,1e-5,1e-3]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=gaussFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = gaussFit(ePoints, *popt1)    # Fit line
# print('Gauss:',popt1)
p00 = np.array([popt1[0],popt1[1],popt1[2],0,0]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=gaussFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = gaussFitPlusBack(ePoints, *popt11)    # Fit line

# print('My G fit:',popt1)
# print('My G+B fit:',popt11)
plt.errorbar(me_a1_4shift_slice['Energy'],me_a1_4shift_slice['Cross-section'],
            yerr=me_a1_4shift_slice['Error'],fmt='.',c='k',label='Data')
plt.plot(ePoints,yfit1,c='r',label='G')
plt.plot(ePoints,yfit11,c='b',label='G+B')
# plt.yscale('log')
plt.text(2.706,0.013,'Pub. Width: 0.1 keV')
# print(type(popt1[1]))
# exit()
w1 = 2.355*popt1[2]*1000
w11 = 2.355*popt11[2]*1000
# print('My G width:',w1)
# print('My G+B width:',w11)
plt.text(2.706,0.012,'G Width: %4.3f keV'%w1)
plt.text(2.706,0.011,'G+B Width: %4.3f keV'%w11)
plt.title('Fit of another a$_{1}$ Resonance in $^{27}$Al Data')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Cross-section (arb. units)')
# plt.legend()
# plt.show()


# Lorentz Fit
p0 = np.array([2.71,1e-2,1e-5]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=lorenFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = lorenFit(ePoints, *popt1)    # Fit line

p00 = np.array([popt1[0],popt1[1],popt1[2],popt11[3],popt11[4]]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=lorenFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = lorenFitPlusBack(ePoints, *popt11)    # Fit line

# print('My L fit:',popt1)
# print('My L+B fit:',popt11)
plt.plot(ePoints,yfit1,c='g',label='L')
plt.plot(ePoints,yfit11,c='m',label='L+B')
w1 = popt1[1]*1000
w11 = popt11[1]*1000
# print('My L width:',w1)
# print('My L+B width:',w11)
plt.text(2.706,0.01,'L Width: %4.3f keV'%w1)
plt.text(2.706,0.009,'L+B Width: %4.3f keV'%w11)
plt.legend()
plt.xlim(2.705,2.72)
# plt.yscale('log')
plt.savefig('resWidthAnalOther.png',dpi=900)
# plt.show()
plt.clf()
# exit()





# Weighted fit of Nelson a1 data
me_a1_4shift_slice = nel_a1_120.query('Energy >= 2.705 and Energy <= 2.72')

# Gauss Fit
ePoints = np.linspace(min(me_a1_4shift_slice['Energy']),max(me_a1_4shift_slice['Energy']),1000)
p0 = np.array([popt1[0],8e-6,2e-3]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=gaussFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = gaussFit(ePoints, *popt1)    # Fit line

p00 = np.array([popt1[0],4e-6,8e-4,0,0]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=gaussFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = gaussFitPlusBack(ePoints, *popt11)    # Fit line

# print('Nelson G fit:',popt1)
# print('Nelson G+B fit:',popt11)
plt.errorbar(me_a1_4shift_slice['Energy'],me_a1_4shift_slice['Cross-section'],
            yerr=me_a1_4shift_slice['Error'],fmt='.',c='k',label='Data')
plt.plot(ePoints,yfit1,c='r',label='G')
plt.plot(ePoints,yfit11,c='b',label='G+B')
# plt.yscale('log')
plt.text(2.706,0.002,'Pub. Width: 0.1 keV')
# print(type(popt1[1]))
# exit()
w1 = 2.355*popt1[2]*1000
w11 = 2.355*popt11[2]*1000
# print('My G width:',w1)
# print('My G+B width:',w11)
plt.text(2.706,0.0018,'G Width: %4.3f keV'%w1)
plt.text(2.706,0.0016,'G+B Width: %4.3f keV'%w11)
plt.title('Fit of another a$_{1}$ Resonance in $^{27}$Al Nelson Data 120$^{o}$')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Cross-section (arb. units)')
# plt.legend()
# plt.show()
# plt.clf()
# exit()


# Lorentz Fit
p0 = np.array([popt11[0],popt11[2],popt11[1]]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=lorenFit, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p0, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit1 = lorenFit(ePoints, *popt1)    # Fit line
print(popt1)
p00 = np.array([popt1[0],popt1[1],popt1[2],popt11[3],popt11[4]]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=lorenFitPlusBack, xdata=me_a1_4shift_slice['Energy'],
    ydata=me_a1_4shift_slice['Cross-section'],p0=p00, sigma=me_a1_4shift_slice['Error'],
    absolute_sigma=True)
yfit11 = lorenFitPlusBack(ePoints, *popt11)    # Fit line
print(popt11)
# print('My L fit:',popt1)
# print('My L+B fit:',popt11)
plt.plot(ePoints,yfit1,c='g',label='L')
plt.plot(ePoints,yfit11,c='m',label='L+B')
w1 = popt1[1]*1000
w11 = popt11[1]*1000
# print('Nelson L width:',w1)
# print('Nelson L+B width:',w11)
plt.text(2.706,0.0014,'L Width: %4.3f keV'%w1)
plt.text(2.706,0.0012,'L+B Width: %4.3f keV'%w11)
plt.xlim(2.705,2.72)
plt.legend()
# plt.yscale('log')
plt.savefig('resWidthAnalNelson120Other.png',dpi=900)
# plt.show()









exit()











# Weighted fit of Nelson a1 data
nel_a1_120 = nel_a1_120.query('Energy >= 2.735 and Energy <= 2.745')

p0 = np.array([2.7405,1e-4,0.001]) # initial guess of parameters
popt1, pcov1 = curve_fit(f=gaussFit, xdata=nel_a1_120['Energy'],
    ydata=nel_a1_120['Cross-section'],p0=p0, sigma=nel_a1_120['Error'],
    absolute_sigma=True)
yfit1 = gaussFit(nel_a1_120['Energy'], *popt1)    # Fit line

p00 = np.array([2.7405,1e-4,0.001,1,1]) # initial guess of parameters
popt11, pcov11 = curve_fit(f=gaussFitPlusBack, xdata=nel_a1_120['Energy'],
    ydata=nel_a1_120['Cross-section'],p0=p00, sigma=nel_a1_120['Error'],
    absolute_sigma=True)
yfit11 = gaussFitPlusBack(nel_a1_120['Energy'], *popt11)    # Fit line

print('Nelson G fit:',popt1)
print('Nelson G+B fit:',popt11)
plt.errorbar(nel_a1_120['Energy'],nel_a1_120['Cross-section'],
            yerr=nel_a1_120['Error'],fmt='.',c='k',label='Nel120 Data')
plt.plot(nel_a1_120['Energy'],yfit1,c='r',label='G')
plt.plot(nel_a1_120['Energy'],yfit11,c='b',label='G+B')
# plt.yscale('log')
plt.title('Fit of an a$_{1}$ Resonance in Nelson 27Al Data 120$^{o}$')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Cross-section (arb. units)')
plt.legend()
plt.show()





exit()
print('Done')


me_p1 = pd.read_table('rMatrix/27Al_rMatrix_p1_angInt.dat',names=colNames)
nel_p1 = pd.read_table('rMatrix/Nelson_pp1_EXFOR.dat',names=colNames)

# Nelson p1 angles: 90, 105, 135, 160
nel_p1_90 = nel_p1.query('Angle == 90')
nel_p1_105 = nel_p1.query('Angle == 105')
nel_p1_135 = nel_p1.query('Angle == 135')
nel_p1_160 = nel_p1.query('Angle == 160')

# plt.scatter(me_p1['Energy'],me_p1['Cross-section'],s=4,c='r',label='Me')
# plt.scatter(nel_p1_90['Energy'],nel_p1_90['Cross-section'],s=4,c='k',label='Me')
# plt.xlim(2.26,2.95)
# plt.yscale('log')
# plt.show()
