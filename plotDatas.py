import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read in the data
colNames = ['Energy','Angle','Cross-section','Error']
me_a1 = pd.read_table('rMatrix/27Al_rMatrix_a1.dat',names=colNames,sep='\s+')
me_p1 = pd.read_table('rMatrix/27Al_rMatrix_p1.dat',names=colNames,sep='\s+')
me_p2 = pd.read_table('rMatrix/27Al_rMatrix_p2.dat',names=colNames,sep='\s+')
nel_a0 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)
nel_a1 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)
nel_p1 = pd.read_table('rMatrix/Nelson_pp1_EXFOR.dat',names=colNames)
nel_p2 = pd.read_table('rMatrix/Nelson_pp2_EXFOR.dat',names=colNames)

# Shift Nelson energies by -4 keV
nel_a0_shift = nel_a0.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_a1_shift = nel_a1.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_p1_shift = nel_p1.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_p2_shift = nel_p2.apply(lambda x: x-.004 if x.name == 'Energy' else x)


me_a1_60 = me_a1.query('Angle==60')
nel_a1_60_shift = nel_a1_shift.query('Angle==120')

plt.errorbar(x=me_a1_60['Energy'],y=me_a1_60['Cross-section'],yerr=me_a1_60['Error'],fmt='ok',alpha=.5,label='me')
plt.errorbar(x=nel_a1_60_shift['Energy'],y=nel_a1_60_shift['Cross-section'],yerr=nel_a1_60_shift['Error'],fmt='ob',alpha=.5,label='nelson')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Differential Cross-section (b/sr)')
plt.title('a1 - 120$^{o}$ Detector')
plt.savefig('plots/a1_120_Detector.png')
plt.legend()
# plt.show()
plt.clf()


me_p1_90 = me_p1.query('Angle==90')
nel_p1_90_shift = nel_p1_shift.query('Angle==90')

plt.errorbar(x=me_p1_90['Energy'],y=me_p1_90['Cross-section'],yerr=me_p1_90['Error'],fmt='ok',alpha=.5,label='me')
plt.errorbar(x=nel_p1_90_shift['Energy'],y=nel_p1_90_shift['Cross-section'],yerr=nel_p1_90_shift['Error'],fmt='ob',alpha=.5,label='nelson')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Differential Cross-section (b/sr)')
plt.title('p1 - 90$^{o}$ Detector')
plt.savefig('plots/p1_90_Detector.png')
plt.legend()
# plt.show()
plt.clf()


me_p1_75 = me_p1.query('Angle==75')
nel_p1_75_shift = nel_p1_shift.query('Angle==105')

plt.errorbar(x=me_p1_75['Energy'],y=me_p1_75['Cross-section'],yerr=me_p1_75['Error'],fmt='ok',alpha=.5,label='me')
plt.errorbar(x=nel_p1_75_shift['Energy'],y=nel_p1_75_shift['Cross-section'],yerr=nel_p1_75_shift['Error'],fmt='ob',alpha=.5,label='nelson')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Differential Cross-section (b/sr)')
plt.title('p1 - 105$^{o}$ Detector')
plt.savefig('plots/p1_105_Detector.png')
plt.legend()
# plt.show()
plt.clf()


me_p2_90 = me_p2.query('Angle==90')
nel_p2_90_shift = nel_p2_shift.query('Angle==90')

plt.errorbar(x=me_p2_90['Energy'],y=me_p2_90['Cross-section'],yerr=me_p2_90['Error'],fmt='ok',alpha=.5,label='me')
plt.errorbar(x=nel_p2_90_shift['Energy'],y=nel_p2_90_shift['Cross-section'],yerr=nel_p2_90_shift['Error'],fmt='ob',alpha=.5,label='nelson')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Differential Cross-section (b/sr)')
plt.title('p2 - 90$^{o}$ Detector')
plt.savefig('plots/p2_90_Detector.png')
plt.legend()
plt.show()
plt.clf()


me_p2_75 = me_p2.query('Angle==75')
nel_p2_75_shift = nel_p2_shift.query('Angle==105')

plt.errorbar(x=me_p2_75['Energy'],y=me_p2_75['Cross-section'],yerr=me_p2_75['Error'],fmt='ok',alpha=.5,label='me')
plt.errorbar(x=nel_p2_75_shift['Energy'],y=nel_p2_75_shift['Cross-section'],yerr=nel_p2_75_shift['Error'],fmt='ob',alpha=.5,label='nelson')
plt.xlabel('E$_{p}$ (MeV)')
plt.ylabel('Differential Cross-section (b/sr)')
plt.title('p2 - 105$^{o}$ Detector')
plt.savefig('plots/p2_105_Detector.png')
plt.legend()
plt.show()
plt.clf()
