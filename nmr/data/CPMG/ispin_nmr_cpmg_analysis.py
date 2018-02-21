'''
File:  ispin_nmr_cpmg_analysis.py

Last update:  9/22/2013: Tycho Sleator -

#############################################################################
Program to plot and analyze CPMG data, in which just a few points (n) are
taken from the top of each echo.  The data are generated in the Spincore
LabView interface from the specification of the number of points per echo, n.
At the beginning of the FID, and in the middle of each echo, n points are
collected at the rate given by the parameter 'SW' in the LabVeiw interface,
(or 'bw' in this program).  The n points from each echo are averaged, and the
average is plotted as a function of the time of the given echo.  The data are
fit to a decaying exponential, from which the value of T2 can be extracted.
#############################################################################
'''

import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.optimize

'''
###########################################################################
Set the folowing parameters based on values from SpinCore LabView interface
###########################################################################
'''
fname = "cpmgTau3750-5pts.csv" # This is the filename containing the data
tau = 3750e-6     # Obtained from value of tau in SpinCore LabView interface
numave = 5            # Number of points/echo, or number of adjacent points to
                      # average.  Obtained from SpinCore LabView interface
T2_0 = 0.01           # Initial guess for T2 (seconds)
'''
###########################################################################
Set the above parameters based on values from SpinCore LabView interface
###########################################################################
'''

'''
The following reads in the data from the specified file (fname)
'''
infile = open(fname,"r")
text = infile.read()      # read file into a string
infile.close()

#index = text.find("@SW=") # Find the position of "@SW="
#text2 = text[index:-1]    # Create new string beginning at "@SW="
#index2 = text2.find('\n') # Find the next CR in this string
# print 'text=',text2[0:index2]
bw = 30000.000 #float(text2[4:index2]) # This gives the bandwidth
# print 'bw = ',bw
# print '1/bw = ',1/bw  # Note that the time interval between points is 1/bw

# Read the data from the the file starting on line 2
s1 = mlab.csv2rec(fname, skiprows=0)
s1 = s1[:-4] #because we have weird size mismatch

t =  (1/bw)*np.arange(len(s1)/2)  #time data

# assign the data to variables with shorter names
s = s1['amplitude__plot_0']
rs = s.reshape(-1,2)
rtp= np.transpose(rs) # rtp[0] is the real part and rtp[1] is the imaginary part of the data

za = rtp[0] + rtp[1]*1j # create complex array
# print 'abs(za)=',abs(za)
# print 'len(za) = ',len(za)

cpmgmag = abs(za)  # get the magnitude of the cpmg data

'''
Average over each set of 'numave' consecutive points in the data.  'numave' is
the number of points per echo (from LabView interface).  This produces data
consisting of Necho points (including the FID).  These data will be fit to an
exponential function to determine the value of T2.
'''

za_reshape = abs(za.reshape(len(za)/numave,numave))
expdata = np.apply_along_axis(sum,1,za_reshape)/numave
print
print 'expdata=',expdata
print 'len(expdata)=', len(expdata)
print
texpdata = tau*np.arange(len(expdata))  # times of each echo

'''
Fit the data to an exponential function: a exp(-t/T2) + c
'''
# Define nonlinear fitting function.  This is the function we fit to the data.
def f(t, A, T2, y0):
    return A*np.exp(-t/T2) + y0

# initial guesses of the parameters: modify these if fit doesn't converge.
A0 = expdata[0];  # Value of initial point (FID amplitude)
T20 = 0.02;
y00 = 0;

# This is the function that does the nonlinear fit:
nlfit, nlpcov = scipy.optimize.curve_fit(f, texpdata, expdata, p0=[A0, T20, y00], sigma=None)

# These are the parameter estimates (assigned to more readible names):
A_nlfit = nlfit[0]
T2_nlfit = nlfit[1]
y0_nlfit = nlfit[2]

'''
Below are the uncertainties in the estimates:  (note that "nlpcov" is the
"covariance matrix".  The diagonal elements of the covariance matrix (nlpcov[i][i])
are the variances of the fit parameters, whose square roots give the standard
deviation or uncertainties in the fit parameters.  The off-diagonal elements give
the correlations between fit parameters.
'''
Asig_nlfit =  np.sqrt(nlpcov[0][0])
T2sig_nlfit = np.sqrt(nlpcov[1][1])
y0sig_nlfit = np.sqrt(nlpcov[2][2])

print ''
print '================='
print ' Fit Parameters  '
print '====================================================='
print 'T2 = ',T2_nlfit,  '+/- ', T2sig_nlfit
print 'A  = ',A_nlfit,    '+/- ', Asig_nlfit
print 'y0 = ',y0_nlfit,  '+/- ', y0sig_nlfit
print '====================================================='
print ''


# create the figure
fig1   = plt.figure(1, figsize = (8,9) )
sp1    = fig1.add_subplot(211)  # this will show that time data
sp2    = fig1.add_subplot(212)  # this will show the frequency data


'''
Plot the data obtained from the iSpin Labview interface
'''
plt.figure(1, figsize = (6,5) )

# draw x and y axes
sp1.axhline(color ='r')
sp1.axvline(color ='r')
#print 'len(az.real)=',len(za.real)
#print 'len(t)=',len(t)

# plot the points
sp1.plot(t,za.real, 'bo-')  # plot the real part (blue)
sp1.plot(t,za.imag, 'ro-')  # plot the imaginary part (red)
sp1.plot(t,abs(za),'ko-')   # plot the absolute value (black)

# label the axes
sp1.set_xlabel('Time (sec)',fontsize=14)
sp1.set_ylabel('iSpin CPMG Data',fontsize=14)

# specify the plot limits
sp1.set_xlim(t[0],t[-1])

# display the figure
#sp1.show()



# plot echo size as a function of time from original pi/2 pulse
#plt.figure(2, figsize = (8,5) )

# draw x and y axes
sp2.axhline(color ='r')
sp2.axvline(color ='r')
#print 'len(az.real)=',len(za.real)
#print 'len(t)=',len(t)

tfit = np.linspace(0,texpdata[-1],num=50)
expfit = f(tfit,A_nlfit,T2_nlfit,y0_nlfit)
sp2.plot(tfit,expfit,'b-')
sp2.plot(texpdata,expdata, 'bo')  # plot the data points (blue)

sp2.text(0.3*texpdata[-1], 0.7*expdata[0],
     '$T_2 \, =\, {0:6.4f} \pm {1:6.5f}$ s '.format(T2_nlfit, T2sig_nlfit), ha='left', va='bottom', size='large')


# label the axes
sp2.set_xlabel('Echo Time (sec)',fontsize=14)
sp2.set_ylabel('Echo Magnitude',fontsize=14)

# specify the plot limits
xrange = texpdata[-1]-texpdata[0]
sp2.set_xlim(texpdata[0]-0.03*xrange,texpdata[-1]+0.03*xrange)
#plt.xlim({texpdata[0],texpdata[-1]}+0.01*xrange*{0,1})
#sp2.set_ylim(0,expdata[0])

# display the figure
plt.show()
