'''
This routine reads in a spincore ".txt" file, and plots the data as well as the
fft of the data in that file

Last update:  2/20/2012, 10/7/2012 by Tycho Sleator
Ammended Feb 2018 by Kaitlyn and Patrick
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
Array y is considered a function of array x
We integrate y as function of x in the range of x1 to x2,
'''
def intsimps(y,x,x1,x2):
    x0 = x[0]
    xmax = x[-1]
    dx = x[1] - x[0]
    i1 = int( np.floor((x1-x0)/dx) ) # index of starting value of x
    i2 = int( np.floor((x2-x0)/dx) ) # index of ending value of x
#    print 'i1 = ',i1, ' : i2 = ',i2
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
    #print 'text=',text2[0:index2]
    bw = float(text2[4:index2]) # This gives the bandwidth
    print 'bw = ',bw
    print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

    # Read the data from the the file starting on line 13
    s1 = mlab.csv2rec(fname, skiprows=12)
    npts = len(s1)/2  # number of complex data points

    print 'npts = ',npts

    t =  (1/bw)*np.arange(npts)  #time data

    # assign the data to variables with shorter names
    s = s1['data']
    rs = s.reshape(-1,2)
    rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary part of the data

    return (t,(rtp[0] + rtp[1]*1j)) # create complex array

# rotate data in complex plane
#theta = -0.0*np.pi
#za = np.exp(1j*theta)*(rtp[0] + rtp[1]*1j) # create complex array


'''
Begin Execution Here
'''

# The following is the filename containing the data
directory = "data/fid/"
filename = "fid0.txt"
fname = directory+filename
print "filename = ",fname
plot_title='Free Induction Decay with\n Pulse Time = 5 $\mu s$\n Spectral Width = 100 kHz\n Spectrometer Frequency = 21.15 MHz'

# read data from file
(t,za) = read_data_file(fname)
bw = 1/(t[1]-t[0])
npts = len(t)

'''
Between the comment signs, we do things to the time signal before taking the fft

npts = 2*npts
za = np.append(za,0*za)
#za = np.append(za,0*za)
t =  (1/bw)*np.arange(npts)  #time data

sigma = 1.8/1000.

za = za*np.exp(-t**2/(2*sigma**2))


Between the comment signs, we do things to the time signal before taking the fft
'''

fza = np.fft.fftshift(np.fft.fft(za))  # take the fft to get the frequency spectrum
f = (bw/npts)*np.arange(npts)-bw/2     # array of frequencies

# Integrate over the frequency data

fintegral = intsimps(fza,f,1709-10000,1709+10000)
print 're(freq integral)/10**12 = ', fintegral.real/10**12
print 'im(freq integral)/10**12 = ', fintegral.imag/10**12

# get maximum value of the data
maxza = np.max([max(za.real),max(za.imag)])  # maximum value

pow = np.floor(np.log10(maxza))  # next power of 10 less than max

'''
CREATE THE FIGURE
'''
# fig1   = plt.figure(figsize=(8,10))
ax1 = plt.axes()

plt.title('{0}'.format(plot_title))

#ax1    = fig1.add_subplot(211)  # this will show that time data
# ax2    = fig1.add_subplot(212)  # this will show the frequency data

'''
Plot the Time Data:
'''
# draw x and y axes
ax1.axhline(color ='k')
ax1.axvline(color ='k')
print 'len(az.real)=',len(za.real)
print 'len(t)=',len(t)

tscale = 1000.0   # change time units to msec
tunits = 'msec'
fscale = 1/tscale
funits = 'kHz'

# plot the points
ax1.plot(t*tscale,za.real/10**pow, '-b', label='Real Part')  # plot the real part (blue)
ax1.plot(t*tscale,za.imag/10**pow, '-r', label='Imaginary Part')  # plot the imaginary part (red)
ax1.plot(t*tscale,np.sqrt(za.real**2 + za.imag**2)/10**pow, '-k', label='Magnitude') # plot the magnitude (black)

# label the axes and display legend
ax1.set_xlabel('Time ('+np.str(tunits)+')',fontsize=14)
ax1.set_ylabel('Signal (x $10^'+str(int(pow))+ '$)',fontsize=14)
ax1.legend()

# specify the plot limits
ax1.set_xlim(t[0]*tscale,t[-1]*tscale)

'''
Plot the Frequency Data:
'''
# draw x and y axes
# ax2.axhline(color ='k')
# ax2.axvline(color ='k')

# # plot the points
# ax2.plot(f*fscale,fza.real, '-b', label='Real Part')  # plot the real part (blue)
# ax2.plot(f*fscale,fza.imag, '-r', label='Imaginary Part')  # plot the imaginary part (red)
# ax2.plot(f*fscale,np.sqrt(fza.real**2 + fza.imag**2), '-k', label='Magnitude')  # plot the magnitude (black)

# # label the axes and display the legend
# ax2.set_xlabel('Frequency ('+np.str(funits)+')',fontsize=14)
# ax2.set_ylabel('Signal',fontsize=14)
# ax2.legend()

# # specify the plot limits
# [fmin,fmax] = [f[0]/20,f[-1]/20]
# ax2.set_xlim(fmin*fscale,fmax*fscale)


'''
Display the Figure
'''
plt.show()
