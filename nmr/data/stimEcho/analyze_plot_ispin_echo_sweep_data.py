'''
File:  'analyze_plot_ispin_echo_sweep_data.py'
This file is derived from the file 'manipulate_plot_ispin_nmr_data.py'.
(which was previously updated:  2/20/2012 by Tycho Sleator)

This file reads in a sequence of spincore ".txt" files generated by the 
"hahn_echo" sequence, extracts numbers from each of those files, generates 
and prints an array from these numbers and plots the result.

The array can be pasted into a fitting program, such as 
'NMR_diffusion_fit_array.py' for further analysis.

Last update:  10/07/2012 at 12:47 pm by Tycho Sleator
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.integrate as integrate

'''
THE FOLLOWING ARE FUNCTION DEFINITIONS
'''
'''
Array y is considered a function of array x
Integrate y as function of x in the range of x1 to x2, 
'''
def intsimps(y,x,x1,x2):
    x0 = x[0]
    xmax = x[-1]
    dx = x[1] - x[0]  
    i1 = int(np.floor((x1-x0)/dx))  # index of starting value of x
    i2 = int(np.floor((x2-x0)/dx)) # index of ending value of x
    print 'i1 = ',i1, ' : i2 = ',i2
    print len(y), len(x)
    print x
    return integrate.simps(y[i1:i2],None,dx)

'''
This reads in complex data from a file and outputs a tuple (t,za),
where t is the time data and za is an array of complex numbers 
corresponding to the time data
'''
def read_data_file(fname):
    infile = open(fname,"r")
    text = infile.read()      # read file into a string
    infile.close()

    index = text.find("@SW=") # Find the position of "@SW="
    text2 = text[index:-1]    # Create new string beginning at "@SW="
    index2 = text2.find('\n') # Find the next CR in this string
    # print 'text=',text2[0:index2]
    bw = float(text2[4:index2]) # This gives the bandwidth
    # print 'bw = ',bw
    # print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

    # Read the data from the the file starting on line 13
    s1 = mlab.csv2rec(fname, skiprows=12)  
    npts = len(s1)/2  # number of complex data points

    # print 'npts = ',npts

    t =  (1/bw)*np.arange(npts)  #time data

    # assign the data to variables with shorter names
    s = s1['data']
    rs = s.reshape(-1,2) 
    rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary 
                          # part of the data

    return (t,(rtp[0] + rtp[1]*1j)) # create complex array

# rotate data in complex plane 
#theta = -0.0*np.pi
#za = np.exp(1j*theta)*(rtp[0] + rtp[1]*1j) # create complex array

'''
EXECUTION BEGINS HERE
'''
# The following is the root filename containing the data
rootname = "stim_echo_sweep"
taumin = 5 # minimum value of tau
Ntaus = 10  # of different values of tau (delays from tau to Ntau*tau)
Delta_tau = 5
taus = taumin + Delta_tau * np.arange(Ntaus) # list of values of tau
#taus = tau*np.arange(1,Ntaus+1)
fnums = range(Ntaus)  # file numbers (range from 0 to 19)

print 'taus = ',list(taus)
print 'fnums = ',fnums

echos = [];  #initialize array of echosizes

for fnum in fnums:    
    fname = rootname+str(100+fnum)[1:]+'.txt'  # generate the filename
    print 'filename = ', fname
	
    # read data from file
    (t,za) = read_data_file(fname)
    bw = 1/(t[1]-t[0])
    npts = len(t)
    
    echosize = intsimps(abs(za),t,0.017,0.023) # compute the 'size' of the echo
    echos.append(echosize)  # add to the array of echo sizes
    print 'echosize =', echosize
rechos = echos/(max(echos))
print 'taus = ', list(taus)
print 'echo sizes = ',list(echos)
print 'relative echo sizes = ',list(rechos)
           
'''
CREATE THE FIGURE
'''
fig1   = plt.figure(figsize=(8,5))
# draw x and y axes
plt.axhline(color ='k')
plt.axvline(color ='k')
#
# plot the points
plt.plot(taus,rechos, 'ob')  # plot the real part (blue)
#
## label the axes
plt.xlabel('Echo Delay $\\tau$ (ms)',fontsize=14)
plt.ylabel('Relative Echo Size' ,fontsize=14)
#
# specify the plot limits
plt.ylim(0,1.05)
#
#Display the Figure
plt.show()                                      
