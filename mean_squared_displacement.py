#############################################
# Plot mean squared displacement using data generated in fortran code
#############################################

def linear(x, a, b):
    return a*x + b

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize

data = np.loadtxt('Mean_squared_displacement_Ar.txt')

DT = 10.0E-15
natoms = 864
num_runs = 50000
all_steps = np.zeros(300)

# Normalize data
for step_loop in range(0, 300):
    all_steps[step_loop] = step_loop*DT*1E12                 # time in ps
    data[step_loop] = data[step_loop]/natoms/num_runs*(1E20) # mean squared displacement in Angstroms^2

# Fit linear section to calculate diffusion coeficient
popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, all_steps[100:300], data[100:300], p0=[8000.0,0])
perr_linear = np.sqrt(np.diag(pcov_linear))
print("slope = %0.2f (+/-) %0.2f" % (popt_linear[0], perr_linear[0]))
print("intercept= %0.2f (+/-) %0.2f" % (popt_linear[1], perr_linear[1]))

diffusion = popt_linear[0]/6.0 # diffusion in A^2/ps
print(diffusion*100.0*100.0*(1E12)/(1E20)) # Difusion in cm^2/s


plt.figure(figsize=(15,8))
plt.rcParams.update({'font.size': 20})
plt.plot(all_steps,data, 'k',label='Data')
plt.plot([0, 3],[popt_linear[1],3.0*popt_linear[0]+popt_linear[1]], '--b',label='Fit')
plt.axis([0, 3, 0, 5])
plt.xlabel('Time (ps)')
plt.ylabel('Mean squared displacement ($\AA ^2$)')
plt.legend(loc='center right')

plt.savefig('Ar_msd.png', bbox_inches='tight')

plt.show()
