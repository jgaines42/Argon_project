#############################################
# Plot mean squared displacement using data generated in fortran code
#############################################

def linear(x, a, b):
    return a*x +b

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize

data = np.loadtxt('Mean_squared_displacement_Ar.txt')

DT = 10E-15
natoms = 864

num_msd = 400
num_runs = 300000-20000-num_msd
print(num_runs)
num_runs1 = data[data.shape[0]-1]
print(num_runs1)
all_steps = np.zeros(num_msd)
data = data[0:data.shape[0]-1]

# Normalize data
for step_loop in range(0, num_msd):
    all_steps[step_loop] = step_loop*DT*1E12                 # time in ps
    data[step_loop] = data[step_loop]/natoms/num_runs*(1E20) # mean squared displacement in Angstroms^2

# Fit linear section to calculate diffusion coeficient
print(all_steps[250])
popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, all_steps[250:350], data[250:350], p0=[8000.0, 0])
perr_linear = np.sqrt(np.diag(pcov_linear))
print("slope = %0.2f (+/-) %0.2f" % (popt_linear[0], perr_linear[0]))
print("intercept= %0.2f (+/-) %0.2f" % (popt_linear[1], perr_linear[1]))

m = popt_linear[0]
b = popt_linear[1]
diffusion = popt_linear[0]/6.0 # diffusion in A^2/ps
print(diffusion*100.0*100.0*(1E12)/(1E20)) # Difusion in cm^2/s


plt.figure(figsize=(15,8))
plt.rcParams.update({'font.size': 20})
plt.plot(all_steps,data, 'k',label='Data')
plt.plot([0, all_steps[num_msd-1]],[b,b+all_steps[num_msd-1]*m], '--b',label='Fit')
plt.axis([0, 4, 0, 7])
plt.xticks(np.arange(0,4.5, 0.5))
plt.yticks(np.arange(0,7.5, 1))
plt.xlabel('Time (ps)')
plt.ylabel('Mean squared displacement ($\AA ^2$)')
plt.legend(loc='center right')
plt.minorticks_on()
plt.savefig('Ar_msd.png', bbox_inches='tight')

plt.show()
