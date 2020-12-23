#############################################
# Plot mean squared displacement using data generated in fortran code
#############################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize


def linear(x, a, b):
    return a * x + b


data = np.loadtxt('Mean_squared_displacement_Ar.txt')

DT = 10.0E-15
natoms = 864

num_msd = 500

num_runs = data[data.shape[0] - 1]
print(num_runs)

all_steps = np.zeros(num_msd)
data = data[0:data.shape[0] - 1].copy()

# Normalize data
for step_loop in range(0, num_msd):
    all_steps[step_loop] = step_loop * DT * 1E12                 # time in ps
    data[step_loop] = data[step_loop] / natoms / num_runs * (1E20)  # mean squared displacement in Angstroms^2

# Fit linear section to calculate diffusion coeficient
fit_start = 300
fit_end = 499
print(all_steps[fit_start])
print(all_steps[fit_end])
popt_linear, pcov_linear = scipy.optimize.curve_fit(linear, all_steps[fit_start:fit_end], data[fit_start:fit_end],
                                                    p0=[8000.0, 0])
perr_linear = np.sqrt(np.diag(pcov_linear))
print("slope = %0.4f (+/-) %0.2f" % (popt_linear[0], perr_linear[0]))
print("intercept= %0.4f (+/-) %0.2f" % (popt_linear[1], perr_linear[1]))

m = popt_linear[0]
b = popt_linear[1]
diffusion = popt_linear[0] / 6.0  # diffusion in A^2/ps
print(diffusion * 100.0 * 100.0 * (1E12) / (1E20))  # Difusion in cm^2/s

# Plot msd
plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})
plt.plot(all_steps, data, 'k', label='Data')
plt.plot([0, all_steps[num_msd - 1]], [b, b + all_steps[num_msd - 1] * m], '--b', label='Fit')
plt.xlabel('Time (ps)')
plt.ylabel('Mean squared displacement ($\AA ^2$)')
plt.legend(loc='center right')
plt.minorticks_on()
plt.savefig('Ar_msd.png', bbox_inches='tight')

plt.show()

m = (data[fit_end] - data[fit_start]) / (all_steps[fit_end] - all_steps[fit_start])
diffusion = m / 6.0  # diffusion in A^2/ps
diffusion = diffusion * 100.0 * 100.0 * (1E12) / (1E20)
print(diffusion)  # Difusion in cm^2/s
