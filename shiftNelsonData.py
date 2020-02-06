import numpy as np
import pandas as pd

# Read in the data
colNames = ['Energy','Angle','Cross-section','Error']
nel_a0 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)
nel_a1 = pd.read_table('rMatrix/Nelson_pa1_EXFOR.dat',names=colNames)
nel_p1 = pd.read_table('rMatrix/Nelson_pp1_EXFOR.dat',names=colNames)
nel_p2 = pd.read_table('rMatrix/Nelson_pp2_EXFOR.dat',names=colNames)

# Shift Nelson energies by 4 keV
nel_a0_shift = nel_a0.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_a1_shift = nel_a1.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_p1_shift = nel_p1.apply(lambda x: x-.004 if x.name == 'Energy' else x)
nel_p2_shift = nel_p2.apply(lambda x: x-.004 if x.name == 'Energy' else x)

# Saving DataFrame to `.dat` file
np.savetxt('rMatrix/Nelson_pa0_EXFOR_shift.dat', nel_a0_shift.values, fmt='%f')
np.savetxt('rMatrix/Nelson_pa1_EXFOR_shift.dat', nel_a1_shift.values, fmt='%f')
np.savetxt('rMatrix/Nelson_pp1_EXFOR_shift.dat', nel_p1_shift.values, fmt='%f')
np.savetxt('rMatrix/Nelson_pp2_EXFOR_shift.dat', nel_p2_shift.values, fmt='%f')
