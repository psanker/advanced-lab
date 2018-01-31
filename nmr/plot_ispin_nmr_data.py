'''
file:  plot_ispin_nmr_data2.py
This routine reads in a spincore ".txt" file, and plots the data in that file
Unlike the file 'plot_ispin_nmr_data.py' it doesn't automatically rescale the
data to fit between 0 and 1.  This has an advantage when plotting several data
sets on the same plot.  
Last update:  9/25/2012 by Tycho Sleator
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# The following is the filename containing the data
#fname = "../Singlepulse_NMR/SPnmr.txt"
fname = "hahn_echo_sweep_205.txt"
print 'filename = ', fname
infile = open(fname,"r")
text = infile.read()      # read file into a string
infile.close()

index = text.find("@SW=") # Find the position of "@SW="
text2 = text[index:-1]    # Create new string beginning at "@SW="
index2 = text2.find('\n') # Find the next CR in this string
#print 'text=',text2[0:index2]
bw = float(text2[4:index2]) # This gives the bandwidth
print 'bw = ',bw
print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

# Read the data from the the file starting on line 13
s1 = mlab.csv2rec(fname, skiprows=12)  

t =  (1/bw)*np.arange(len(s1)/2)  #time data

# assign the data to variables with shorter names
s = s1['data']
rs = s.reshape(-1,2) 
rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary part of the data

za = rtp[0] + rtp[1]*1j # create complex array


# get maximum value of the data
maxza = np.max([max(za.real),max(za.imag)])  # maximum value

# create the figure
plt.figure(1, figsize = (8,5) )               

# draw x and y axes
plt.axhline(color ='r')
plt.axvline(color ='r')
print 'len(az.real)=',len(za.real)
print 'len(t)=',len(t)


# plot the points
plt.plot(t,np.abs(za),'-k',linewidth=2) # plot the absolute value (black)
plt.plot(t,za.real, '-b')  # plot the real part (blue)
plt.plot(t,za.imag, '-r')  # plot the imaginary part (red)

# label the axes
plt.xlabel('Time (sec)',fontsize=14)
plt.ylabel('Signal',fontsize=14)

# specify the plot limits
plt.xlim(t[0],t[-1])

# display the figure
plt.show()                                      
