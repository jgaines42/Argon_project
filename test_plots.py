def linear(x, a, k):
    return a*x + k

import math # type: ignore
import numpy as np  # type: ignore
import matplotlib.pyplot as plt   # type: ignore
import scipy as scipy
from scipy import optimize
###############################
NVE_start = 10000

dt_NVT = (2.5E-15)*1.0E12          # Timestep in NVT in ps
dt_NVE = (10.0E-15)*1.0E12         # Timestep in NVE in ps
data1 = np.loadtxt('Temp_Ar.txt')  # NVT Temp

data2 = np.arange(data1.shape[0])
data2 = data2*dt_NVE                # Time in NVT

# Plot NVE temperature
plt.figure(figsize=(15,8))
plt.plot(data2,data1, 'k')
plt.xlabel('time (ps)')
plt.ylabel('Temperature (K)')
#plt.axis([0, 100, 86, 100])
plt.plot([0, 100], [np.mean(data1[NVE_start+1:100000]),np.mean(data1[NVE_start+1:100000])], 'r')
plt.plot([0, 100], [94.4, 94.4], '--r')
plt.show()
print(np.mean(data1[0:NVE_start]))
print(np.mean(data1[NVE_start+1:100000]))

# Plot 100 time step moving average of Temp

temp_binned = np.zeros(data1.shape[0])
for i in range(100,data1.shape[0]):
    temp_binned[i] = np.mean(data1[i-100:i])

plt.figure(figsize=(15,8))
plt.plot(data2[100:data1.shape[0]-1], temp_binned[100:data1.shape[0]-1])
print(np.mean(temp_binned[100:data1.shape[0]-1]))
#plt.axis([0, 500, 86, 100])
plt.xlabel('time (ps)')
plt.ylabel('Temperature (K)')
plt.show()

temp_binned = np.zeros(int(data1.shape[0]/100))
for i in range(1,temp_binned.shape[0]):
    temp_binned[i] = np.mean(data1[(i-1)*100+1:(i*100)])
temp_binned[0] = temp_binned[1]
plt.plot(temp_binned)
print(np.mean(temp_binned))
#plt.axis([0, 500, 86, 100])
plt.xlabel('time (ps)')
plt.ylabel('Temperature (K)')
plt.show()
print(temp_binned.shape[0])
np.mean(data1[i-100:i])
data2 = np.arange(temp_binned.shape[0])
popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, data2[500:1000],temp_binned[500:1000], p0=[1,0])
perr_linear = np.sqrt(np.diag(pcov_linear))
print("slope = %0.2f (+/-) %0.2f" % (popt_linear[0], perr_linear[0]))
print("intercept= %0.2f (+/-) %0.2f" % (popt_linear[1], perr_linear[1]))