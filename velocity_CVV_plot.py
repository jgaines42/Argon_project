###################################################
# Plot CVV autocorrelation
############################################

# Define to fit exponential
def exponential(x, a, k):
    return a*np.exp(x*k)

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
data = np.loadtxt('Velocity_Ar.txt')
import scipy as scipy
from scipy import optimize  # type: ignore
import math

average_auto = data
DT = 10.0E-15/(1E-12) #ps
DT_fs = 10.0E-15 #s
NUMBER_ATOMS = 864
NUMBER_TIME = 50000
EQUILIBRIUM_START = 0
NUMBER_STARTS = NUMBER_TIME - EQUILIBRIUM_START
all_steps = np.zeros([len(data)])
average_auto = np.divide(average_auto, NUMBER_STARTS)
average_auto = np.divide(average_auto, NUMBER_ATOMS)

for i in range(0,len(average_auto)):
    all_steps[i] = i*DT

print(average_auto[0])
plt.figure(figsize=(15,8))
plt.rcParams.update({'font.size': 20})
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.plot(all_steps,average_auto, 'k')
plt.plot([0, 1.5], [0, 0])
plt.axis([0, 1.5, -10000,60000])

#plt.show()
plt.savefig('Ar_CVV.png', bbox_inches='tight')

plt.figure(figsize=(15,8))
plt.plot(all_steps,average_auto/average_auto[0], 'k')

plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.plot([0, 1.5], [0, 0])
plt.axis([0, 1.5, -0.2,1])

plt.savefig('Ar_CVV_normalized.png', bbox_inches='tight')
plt.show()


D = np.sum(average_auto[0:150])/3.0 *DT_fs # m^2/s
D = D*(100*100)
print(D)
np.savetxt('Ar_CVV_data.txt', average_auto)

popt_exponential, pcov_exponential = scipy.optimize.curve_fit(exponential, all_steps[50:150], average_auto[50:150], p0=[-8000.0,-2])
perr_exponential = np.sqrt(np.diag(pcov_exponential))
print("pre-exponential factor = %0.2f (+/-) %0.2f" % (popt_exponential[0], perr_exponential[0]))
print("rate constant = %0.2f (+/-) %0.2f" % (popt_exponential[1], perr_exponential[1]))

plt.figure(figsize=(15,8))
plt.plot(all_steps[50:150], exponential(all_steps[50:150], *popt_exponential), 'k--')
plt.plot(all_steps[50:150],average_auto[50:150])
plt.show()

# intergral -> 1/k *exp(kx)
k = popt_exponential[1] # 1/ps
a = popt_exponential[0] # m^2

# Integrate from all_steps[150] to infinity
k = k*1.0E12            # 1/s

integral_15 = np.sum(average_auto[0:150])*DT_fs #m^2/s
integral_15p = 0.0-(1.0/k)*a*math.exp(k*all_steps[149]*1.0E-12)    #m^2/s


D = (integral_15+integral_15p)/3.0 # m^2/s
D = (D)*(100*100)                  # cm^2/s
print(D)