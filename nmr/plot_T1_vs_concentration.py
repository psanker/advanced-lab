'''
Shows the dependence of T1 on sample concentration by plotting the inverse of T1 versus concentration
'''
import numpy as np
import sympy as sm
import matplotlib.pyplot as plt

sm.init_printing(use_latex=True, use_unicode=True)

T_ones = np.array([3846.361, 456.563, 450.840, 141.635, 70.924])
Conc   = np.array([0.0, 1.25, 2.5, 5.0, 10.0])
plt.figure()
plt.plot(Conc, 1./T_ones,'bo',ls='--')
plt.ylabel('1 / $T_1$ (ms$^{-1}$)')
plt.xlabel('Concentration (mM)')
plt.show()
