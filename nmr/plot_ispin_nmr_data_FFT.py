'''
file:  plot_ispin_nmr_data_FFT.py
This routine reads in a spincore ".txt" file, and plots the time data
in that file as well as the frequency data (FFT).  Since the interesting
part of the data is at low frequencies, the frequency range of the plot is
reduced by a scale factor 'sf'.

Last update:  9/25/2012, 10/07/2012 by Tycho Sleator
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# The following is the filename containing the data
fname = "data/fid/fid0.txt"

# The factor by which the frequency scale is expanded around f=0
sf = 1.  # PK CHANGE: Converted sf from int to float

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

t =  (1/bw)*np.arange(len(s1)/2)  # time data
npts = len(t)  # number of data points

# assign the data to variables with shorter names
s = s1['data']
rs = s.reshape(-1,2)
rtp= np.transpose(rs) # rtp[0] is the real part and
                      # rtp[1] is the imaginary part of the data

za = rtp[0] + rtp[1]*1j # create complex array

# get maximum value of the data
# maxza = np.max([max(za.real),max(za.imag)])  # maximum value

# create the figure
fig1   = plt.figure(1, figsize = (8,10) )
ax1    = fig1.add_subplot(211)  # this will show that time data
ax2    = fig1.add_subplot(212)  # this will show the frequency data
'''
Plot the Time Data:
'''
# draw x and y axes
ax1.axhline(color ='k')
ax1.axvline(color ='k')
# print 'len(az.real)=',len(za.real)
# print 'len(t)=',len(t)

tscale = 1000.0   # change time units to msec
tunits = 'msec'
fscale = 1/tscale
funits = 'khz'


# plot the points
ax1.plot(t*tscale,np.abs(za),'-k',linewidth=0.5) # plot the absolute value (black)
ax1.plot(t*tscale,za.real, '-b')  # plot the real part (blue)
ax1.plot(t*tscale,za.imag, '-r')  # plot the imaginary part (red)

# label the axes
ax1.set_xlabel('Time ('+np.str(tunits)+')',fontsize=14)
ax1.set_ylabel('Signal',fontsize=14)

# specify the plot limits
ax1.set_xlim(t[0]*tscale,t[-1]*tscale)

'''
Plot the Frequency Data:
'''
fza = np.fft.fftshift(np.fft.fft(za))  # take fft to get frequency spectrum
f = (bw/npts)*np.arange(npts)-bw/2     # array of frequencies


# draw x and y axes
ax2.axhline(color ='k')
ax2.axvline(color ='k')

# plot the points
ax2.plot(f*fscale,fza.real, '-b')  # plot the real part (blue)
ax2.plot(f*fscale,fza.imag, '-r')  # plot the imaginary part (red)
# plot the magnitude (black)
ax2.plot(f*fscale,np.sqrt(fza.real**2 + fza.imag**2), '-k')

# label the axes
ax2.set_xlabel('Frequency ('+np.str(funits)+')',fontsize=14)
ax2.set_ylabel('Signal',fontsize=14)

# specify the plot limits
[fmin, fmax] = [f[0]/sf, f[-1]/sf]
ax2.set_xlim(fmin*fscale,fmax*fscale)

'''
Display or save the Figure
'''
plt.show()
#plt.savefig(fname.replace(".txt",".png"),format='png')
