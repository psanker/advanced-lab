'''
This file is derived from the file 'plot_ispin_nmr_csv.py'.

This routine reads in the pulsewidth csv file, plots the data, and identifies the location of the pulses.

Amended Feb 2018 by Kaitlyn and Patrick
'''

import numpy as np
import scipy as sp
import sympy as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.integrate as integrate


# Allows LaTeX output in Jupyter and in matplotlib
sm.init_printing(use_latex=True, use_unicode=True)

'''
This reads in complex data from a file and outputs a tuple (t,za,s), where t is the time data, za is an array of complex numbers corresponding to the time data, and s is the amplitude of the raw signal
'''
def read_data_file(fname):
    infile = open(fname,"r")
    text = infile.read()      # read file into a string
    infile.close()

    bw = 100000.000# This gives the bandwidth
    print 'bw = ',bw
    print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

    # Read the data from the the file starting on line 1
    s1 = mlab.csv2rec(fname, skiprows=0)
    npts = len(s1)/2  # number of complex data points

    print 'npts = ',npts

    t =  s1['pulse_width_us__plot_0']  #time data

    # assign the data to variables with shorter names
    s = s1['amplitude__plot_0']
    if len(t) % 2 == 1:
        s = s[1:]
        t = t[1:]
    rs = s.reshape(-1,2)
    rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary part of the data

    return (t,(rtp[0] + rtp[1]*1j),s) # create complex array

# rotate data in complex plane
#theta = -0.0*np.pi
#za = np.exp(1j*theta)*(rtp[0] + rtp[1]*1j) # create complex array


'''
Begin Execution Here
'''

# The following is the filename containing the data
directory = "data/pulseWidth/"
filename = "pulsewidthfind.csv"
fname = directory+filename
print "filename = ",fname
plot_title='Pulse Width'

# read data from file
(t,za,A) = read_data_file(fname)
bw = 1/(t[1]-t[0])
npts = len(t)


# get maximum value of the data
maxza = np.max([max(za.real),max(za.imag)])  # maximum value

pow = np.floor(np.log10(maxza))  # next power of 10 less than max

'''
CREATE THE FIGURE
'''
ax1 = plt.axes()

plt.title('{0}'.format(plot_title))

'''
Plot the Time Data:
'''
# draw x and y axes
ax1.axhline(color ='k')
ax1.axvline(color ='k')
print 'len(az.real)=',len(za.real)
print 'len(t)=',len(t)

tscale = 1e3   # change time units to msec
tunits = 'msec'
pi2_pulse = t[A.argmax()]*tscale
pi_pulse  = 11.25*tscale
pulse_uncertainty = 500*tscale / 2 #500microsec from spincore graph precision

# plot the points
ax1.plot(t*tscale,A,label='Signal')
ax1.axvline(pi2_pulse, ls='--',color='g',alpha=.5,label='$\\pi/2$ pulse = {} $\pm$ {} ms'.format(pi2_pulse, pulse_uncertainty))
ax1.axvline(pi_pulse, ls='--',color='r',alpha=.5,label='$\\pi$ pulse = {} $\pm$ {} ms'.format(pi_pulse, pulse_uncertainty))

# label the axes and display legend
ax1.set_xlabel('Time ('+np.str(tunits)+')',fontsize=14)
ax1.set_ylabel('Signal (x $10^'+str(int(pow))+ '$)',fontsize=14)
ax1.legend(loc='lower right')

# specify the plot limits
ax1.set_xlim(t[0]*tscale,t[-1]*tscale)


'''
Display the Figure
'''
plt.show()
