import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Resonance energies in 'nelson' and my own data 'seb'. Trying to figure out
# relative energy shifts or a calibration I can apply

# Excitation energies in MeV
# Nelson values are quoted from his table of 27Al(p,a1) resonances found in:
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.30.755
'''
`x` means not seen, `-` means seen in my data
2.4418--
2.4716--
2.4883--
2.5172--
2.5728--
2.6044--
2.7132-
2.7410-
2.8099-
2.8226x
2.8458-
2.8765-
2.9110-
2.9207-
3.0171-
3.0397-
'''
# a1 resonances
nelsonLAB = np.array([2.4418,2.4716,2.4883,2.5172,2.5728,2.6044,2.7132,2.7410,2.8099,2.8458,2.8765,2.9110,2.9207,3.0171,3.0397])
sebLAB =    np.array([2.4392,2.4714,2.4842,2.5130,2.5683,2.6010,2.7091,2.7367,2.8062,2.8412,2.8723,2.9065,2.9155,3.0127,3.0345])    # Estimated by eye

e_knots = np.linspace(2.1,3.21,1000)


order=1
z1 = np.polyfit(sebLAB,nelsonLAB,order,full=True)
# print('Residuals Order %s:'%order)
# print(z1[1],'\n\n')
f1 = np.poly1d(z1[0])  # function takes in old sebLAB points and calculates new ones

print('\n\nPolyfit coefficients `Linear`:')
for ind in range(len(z1[0])):
    print('a%d :'%ind,'{:5e}'.format(np.flip(z1[0])[ind]))


order=2
z2 = np.polyfit(sebLAB,nelsonLAB,order,full=True)
# print('Residuals Order %s:'%order)
# print(z2[1],'\n\n')
f2 = np.poly1d(z2[0])  # function takes in old sebLAB points and calculates new ones

print('\n\nPolyfit coefficients `Quadratic`:')
for ind in range(len(z2[0])):
    print('a%d :'%ind,'{:5e}'.format(np.flip(z2[0])[ind]))


order=3
z3 = np.polyfit(sebLAB,nelsonLAB,order,full=True)
# print('Residuals Order %s:'%order)
# print(z3[1],'\n\n')
f3 = np.poly1d(z3[0])  # function takes in old sebLAB points and calculates new ones

print('\n\nPolyfit coefficients `Cubic`:')
for ind in range(len(z3[0])):
    print('a%d :'%ind,'{:5e}'.format(np.flip(z3[0])[ind]))

# plt.scatter(sebLAB,nelsonLAB,s=16,label='Original')
# plt.plot(e_knots,f1(e_knots),label='Linear')
# plt.plot(e_knots,f2(e_knots),label='Quad')
# plt.plot(e_knots,f3(e_knots),label='Cubic')
# plt.xlabel('Seb')
# plt.ylabel('Nelson')
# plt.title('$^{27}$Al(p,a$_{1}$) - Peak calibration')
# plt.legend()
# plt.show()
# plt.clf()

# # Lab energy values are converted to Excitation energies (keV precision)
# nelsonEX = np.array([13.938,13.967,13.983,14.011,14.064,14.095,14.200])
# sebEX = np.array([13.936,13.966,13.979,14.007,14.061,14.091])    # Estimated by eye




# Read in RMatrix outputs for 27Al data and try recalibrating
colNames = ['Energy','Angle','Cross-section','Error']
channels = ['a1','p1','p2']
for ch in channels:
    df = pd.read_table('rMatrix/27Al_rMatrix_%s_angInt.dat'%ch,names=colNames)

    # Create DataFrames containing new energy calibration for each order
    df1 = df.apply(lambda x: f1(x) if x.name == 'Energy' else x)
    df2 = df.apply(lambda x: f2(x) if x.name == 'Energy' else x)
    df3 = df.apply(lambda x: f3(x) if x.name == 'Energy' else x)

    # Saving DataFrame to `.dat` file
    np.savetxt('rMatrix/27Al_rMatrix_%s_angInt_lin.dat'%ch, df1.values, fmt='%f')
    np.savetxt('rMatrix/27Al_rMatrix_%s_angInt_quad.dat'%ch, df2.values, fmt='%f')
    np.savetxt('rMatrix/27Al_rMatrix_%s_angInt_cubic.dat'%ch, df3.values, fmt='%f')
