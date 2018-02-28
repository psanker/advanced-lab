'''
ispin_nmr_T1_analysis.py:
    
Last updates:  2/21/2012, 10/6/2012, 10/25/2012 - 12:04am by Tycho Sleator

Before, tau_max needed to be computed by "hand" from values of tau_min, 
delta_tau, and n_tau, obtained from the Batch file. Now this is calculated 
automatically.  Improved the output formatting. 

--------------------------------------------------------------------------------
This program determines T1 from a set of spincore ".txt" files, generated by
running the spincore batch file called 'T1_IR_sweep.bat'.  Each ".txt" file 
contains the data from an "inversion recovery" sequence, where a pi pulse 
inverts the signal, and after a variable delay "tau", a pi/2 pulse is applied 
and the resulting free-induction decay (fid) is recorded.  The data in the files 
differ only in the value of tau.  Fitting the size of the fid as a function of 
tau to an exponential decay allows one to extract the spin-lattice relaxation 
time T1.

More specifically, this program:
- plots the data from each of the ".txt" data files.
- Integrates over the resonance line in the frequency spectrum to get the 'size'
  of the signal for each value of tau.
- The integral of the data from each file is rotated in the complex plane to
  make it real, in such a way to compensate for slow phase rotations over time.
- Fits the signal size as a function of tau to a decaying exponental and prints
  the fit parameters (including value of T1), and plots both the data and fit. 

======================================
 Instructions for using this program:
======================================
User needs to set values of:

tau_min     # minimum value of tau (in ms)
delta_tau   # difference between consecutive values of tau (in ms)
n_tau       # number of different values of tau:  this is the number of 
            # files in the directory
These values can be obtained from the corresponding quantities: TAU_MIN, 
DELTA_TAU, and N_TAU in the 'T1_IR_sweep.bat' file.  Note that the units of 
these quantities in the .bat file are microseconds.
    
The range of integration in frequency space (in Hz) for computation of the 
"size" of the FID.  This can be determined by looking at a plot of the FID 
data and the FFT.  

f_min       # minum value of frequency range (in Hz)
f_max       # maximum value of frequency range
    
Initial guesses for the values of the exponential fit perameters that model
the data:  y = A exp(-t/T1) + y0

A0          # initial guess for A
T10         # initial guess for T1
y00         # initial guess for y0

To set these parameters, scroll down to where it says:

    'BEGIN EXECUTION HERE'

shortly below the 'BEGIN EXECUTION HERE', you will see places to enter the 
parameter values.
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.integrate as integrate
import scipy.optimize

import glob

'''
Define exponential function for fit to data
'''
def fitfunc(t, A, T1, y0):         
    return A*np.exp(-t/T1) + y0
    
'''
Array y is considered a function of array x
We integrate y as function of x in the range of x1 to x2, 
'''
def intsimps(y,x,x1,x2):
    x0 = x[0]
    xmax = x[-1]
    dx = x[1] - x[0]  
    i1 = int(np.floor((x1-x0)/dx))  # index of starting value of x
    i2 = int(np.floor((x2-x0)/dx)) # index of ending value of x
    return integrate.simps(y[i1:i2],None,dx)

'''
This reads in complex data from a Spincore ".txt" file and outputs a 
tuple (t,za), where t is the time data and za is an array of complex 
numbers corresponding to the time data
'''
def read_data_file(fname):
    infile = open(fname,"r")
    text = infile.read()       # read file into a string
    infile.close()

    index = text.find("@SW=")  # Find the position of "@SW="
    text2 = text[index:-1]     # Create new string beginning at "@SW="
    index2 = text2.find('\n')  # Find the next CR in this string
    #print 'text=',text2[0:index2]
    bw = float(text2[4:index2]) # This gives the bandwidth

    # Read the data from the the file starting on line 13
    s1 = mlab.csv2rec(fname, skiprows=12)  
    npts = len(s1)/2  # number of complex data points

    t =  (1/bw)*np.arange(npts)  #time data
    s = s1['data']  # signal data
    rs = s.reshape(-1,2) 
    rtp= np.transpose(rs) # rtp[0] is the real part of the data and  
                          # rtp[1] is the imaginary part of the data

    return (t,(rtp[0] + rtp[1]*1j)) # create complex array

'''
This plots the FID data saved by the Spincore software
'''
def plotfid(t, za, f, fza, sf, fname):
    '''
    Create the Figure
    '''
    fig1   = plt.figure(figsize=(8,10))
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
    ax1.plot(t*tscale,za.real, '-b')  # plot the real part (blue)
    ax1.plot(t*tscale,za.imag, '-r')  # plot the imaginary part (red)

    # label the axes
    ax1.set_xlabel('Time ('+np.str(tunits)+')',fontsize=14)
    ax1.set_ylabel('Signal',fontsize=14)

    # specify the plot limits
    ax1.set_xlim(t[0]*tscale,t[-1]*tscale)
    # ax1.set_xlim(0,4)
    '''
    Plot the Frequency Data:
    '''
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
    [fmin,fmax] = [f[0]/sf,f[-1]/sf]
    ax2.set_xlim(fmin*fscale,fmax*fscale)
    '''
    Display or save the Figure
    '''
#    plt.show()   
    # plt.savefig(fname.replace(".txt",".png"),format='png')
    plt.close()  # remove the current figure from memory

'''
Given a complex number, this returns the angle that rotates that complex
number into a real number without changing the sign of the real part. 
'''
def angleRotateReal(z):
    # Given a complex number, this returns the angle that rotates that complex
    # number into a real number, without changing the sign of the real part 
    lr = (1 - 1*np.real(np.sign(z)))/2.0 # 0 if right half-plane and 
                                         # 1 if left half-plane
    angle = np.angle(z)
    if (lr * np.pi - angle) < np.pi:
        return (lr * np.pi - angle)
    else:                             # this statement is not necessary, but it
        return (lr-2) * np.pi - angle # keeps the angles between -180 and +180  
                                       

'''
This takes an array of complex numbers and returns an array of the same size
of real numbers, which can be used for fitting to an exponential.
'''
def complexArrayToReal(zsignal):
    # This creates a real valued array from a complex valued array.
    realsignal = np.array([])

    # print 'zsignal = ', zsignal
    # rotangle is the angle for rotation of first point of zsignal onto negative
    # real axis.
    print('Rotation angles of consecutive data points:')
    print('-------------------------------------------')
    rotangle = np.pi - np.angle(zsignal[0]) # rotate to negative real axis
    if rotangle > np.pi:              # this statement is not necessary, but it
        rotangle = rotangle - 2*np.pi # keeps the angles between -180 and +180
                                      
    print('rotation angle = {: 2.2f} degrees'.format((180/np.pi) * rotangle))

    # rotate all points in the array by this angle.
    zsignal = np.exp(1j * rotangle)*zsignal
    realsignal = np.append(realsignal,zsignal[0].real)
    
    for i in range(1,len(zsignal)): # start with 2nd point (1st already rotated)
    
        # return smallest angle that rotates i'th element to real axis
        rotangle = angleRotateReal(zsignal[i]) 
        print('rotation angle = {: 2.2f} degrees'.format((180/np.pi) * rotangle))

        ## limit rotation angle to 45 degrees (pi/4)        
        #if np.abs(rotangle) > np.pi/4:
        #    rotangle = np.pi/4 * rotangle/abs(rotangle)
        #    print '* corrected angle = ',(180/np.pi) * rotangle, 'degrees'

        # rotate all elements in the array by this amount 
        # (ith element will become real).
        zsignal = np.exp(1j * rotangle) * zsignal 
        
        # append ith element of zsignal to realsignal
        realsignal = np.append(realsignal,zsignal[i].real)
        # print 'realsignal = ', realsignal
    return realsignal


'''
##########################
#  BEGIN EXECUTION HERE  #
##########################
'''
# Put here the values of the tau's from the information in the batch file:
tau_min = 1    # minimum value of tau (in ms)
delta_tau = 200  # difference between consecutive tau values
n_tau = 20      # number of different values of tau:  this is the number of 
                # files in the directory
# tau_max is computed automatically

# The range of integration in frequency space (in Hz) for computation of the 
# "size" of the FID.  
f_min = -400 #(Hz)
f_max = 400

# The factor by which the frequency scale is expanded around f=0 in the 
# frequency plot.
sf = 50

# Initial guess for T1 fit perameters
A0 = -20
T10 = 100 # (ms)
y00 = 10

# The following is the filename containing the data
rootname = "T1_sweep_00"

signal = np.array([])  # array where fid sizes are stored 
zsignal = np.array([]) # complex array
fnumbers = range(0, n_tau)

# The following reads in FID data, takes FFT, plots data and FFT, 
# integrates FFT and appends complex valued integral to array zsignal

print('\nFile Name                Frequency Integral/10^12')
print('---------                ------------------------')

# PK: New directory iteration
file_count = 0

for fname in sorted(glob.glob('*.txt')):
    if file_count >= n_tau:
        continue

    file_count += 1

    # read data from file 
    (t,za) = read_data_file(fname)
    bw = 1/(t[1]-t[0])
    npts = len(t)    

    fza = np.fft.fftshift(np.fft.fft(za))  # take fft to get frequency spectrum
    f = (bw/npts)*np.arange(npts)-bw/2     # array of frequencies

    # display a plot of the time and frequency data, or save plot to a pdf file
    # t is the time axis, za is the signal as a function of time
    # f is the frequency axis, fza is the fft of za
    # sf is the scale factor for the frequency, fname is the file name.  
    plotfid(t, za, f, fza, sf, fname) 
    
    # Integrate over the frequency data
    fintegral = intsimps(fza,f,f_min,f_max)
    print(fname,': {: 1.3e} + {: 1.3e}j'.format(fintegral.real/10**12,fintegral.imag/10**12))
        
    # Add the integrated signal to the "signal" vs "tau" data
    signal = np.append(signal, fintegral.real/10**12)
    zsignal = np.append(zsignal, fintegral/10**12)
    
# OLD
#for fnumber in fnumbers:
#    fname = rootname+str(100 + fnumber)[1:]+".txt"
#
#    # read data from file 
#    (t,za) = read_data_file(fname)
#    bw = 1/(t[1]-t[0])
#    npts = len(t)    
#
#    fza = np.fft.fftshift(np.fft.fft(za))  # take fft to get frequency spectrum
#    f = (bw/npts)*np.arange(npts)-bw/2     # array of frequencies
#
#    # display a plot of the time and frequency data, or save plot to a pdf file
#    # t is the time axis, za is the signal as a function of time
#    # f is the frequency axis, fza is the fft of za
#    # sf is the scale factor for the frequency, fname is the file name.  
#    plotfid(t, za, f, fza, sf, fname) 
#    
#    # Integrate over the frequency data
#    fintegral = intsimps(fza,f,f_min,f_max)
#    print fname,': {: 1.3e} + {: 1.3e}j'\
#        .format(fintegral.real/10**12,fintegral.imag/10**12)
#        
#    # Add the integrated signal to the "signal" vs "tau" data
#    signal = np.append(signal, fintegral.real/10**12)
#    zsignal = np.append(zsignal, fintegral/10**12)

print("\n complex t1 data:\n-------------------\n ",zsignal)
print("number of points = ",len(signal))
print("\n")

realsignal = complexArrayToReal(zsignal)

print("\n real t1 data:\n-------------------\n ",realsignal)

# values of tao for which data were taken
tau_max = tau_min + delta_tau * (n_tau - 1)  # maximum value of tau 
tau = np.linspace(tau_min, tau_max, n_tau) 

# use this array of tau's to plot the fit
taufit = np.linspace(0, 1.05*tau_max, 5*n_tau)

# This is the function that does the nonlinear fit:
nlfit, nlpcov = \
    scipy.optimize.curve_fit(fitfunc, tau, realsignal, p0=[A0, T10, y00])

# Best values of fit parameters 
A_nlfit = nlfit[0]
T1_nlfit = nlfit[1]
y0_nlfit = nlfit[2]

# Uncertainties in the parameters
Asig_nlfit =  np.sqrt(nlpcov[0][0])
T1sig_nlfit = np.sqrt(nlpcov[1][1])
y0sig_nlfit = np.sqrt(nlpcov[2][2])

print('\n')
print('=================================')
print('       Fit Paremeters    ')
print('  T1 = {: 1.2f}   +- {:1.2f} ms'.format(T1_nlfit,T1sig_nlfit))
print('  A  = {: 1.2e} +- {:1.1e}'.format(A_nlfit,Asig_nlfit))
print('  y0 = {: 1.2e} +- {:1.1e}'.format(y0_nlfit,y0sig_nlfit))
print('=================================')

# draw the figure of signal vs. tau
final = plt.figure(figsize=(8,5))
plt.axhline(color ='k')  # draw x and y axes
plt.axvline(color ='k')
plt.plot(np.linspace(tau_min,tau_max,n_tau),realsignal, 'ob') # plot the pts
plt.plot(taufit,fitfunc(taufit,A_nlfit,T1_nlfit,y0_nlfit),'-r') # plot the fit
plt.xlabel('Delay Time, $\\tau$ (ms)',fontsize=14)  # x-axis label
plt.ylabel('FID Amplitude',fontsize=14)    # y-axis label
# print value of T1
plt.text(0.7*tau_max, 0.1*max(realsignal) + 0.9*min(realsignal),
'$y \,=\, y_0+ A\, \exp(-t/T_1)$\n\n$T_1 \,=\, {:.3f}$ ms\n $A \,=\, {:.3f}$ \
    \n $y_0 \,=\, {:.3f}$'.format(T1_nlfit, A_nlfit,y0_nlfit),fontsize=14)

# PK: change final -> plt
plt.show()

# OLD
# final.show()
