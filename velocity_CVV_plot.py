###################################################
# Plot CVV autocorrelation using data generated in fortran code
############################################
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import scipy as scipy
from scipy import optimize  # type: ignore
import math


# Define to fit exponential
def exponential(x, a, k):
    return a * np.exp(x * k)


# Load file
# Velocity_Ar.txt is the autocorrelation data from the run
data = np.loadtxt('Velocity_Ar.txt')

DT = 10.0E-15 / (1E-12)       # Time step in ps
DT_fs = 10.0E-15            # Time step in s
NUMBER_ATOMS = 864          # Number of atoms
NUMBER_TIME = 400000         # Number of time steps
EQUILIBRIUM_START = 0       # When equlibirum starts (in the NUMBER_TIME)

NUMBER_STARTS = data[data.shape[0] - 1]
print(NUMBER_STARTS)
data = data[0:data.shape[0] - 1].copy()
print(data.shape[0])
print(len(data))
average_auto = data

# Normalilze data by the number of time steps in it and the number of atoms
average_auto = np.divide(average_auto, NUMBER_STARTS)
average_auto = np.divide(average_auto, NUMBER_ATOMS)

# Make time data
all_steps = np.zeros([len(data)])
for i in range(0, len(average_auto)):
    all_steps[i] = i * DT         # in ps

# Plot CVV(t)
plt.figure(figsize=(15, 8))
plt.rcParams.update({'font.size': 20})
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.plot(all_steps, average_auto, 'k')
plt.plot([0, 3], [0, 0])
plt.xticks(np.arange(0, 3.5, 0.5))

plt.axis([0, 3, -10000, 60000])


plt.savefig('Ar_CVV.png', bbox_inches='tight')
plt.show()
print(average_auto[1])

# Plot normalized by dt = 0
plt.figure(figsize=(15, 8))
plt.plot(all_steps, average_auto / average_auto[0], 'k')
plt.xlabel('time (ps)')
plt.ylabel('$C_{vv}(t), (\mathrm{m^2/s^2})$')
plt.plot([0, 3], [0, 0])
plt.axis([0, 3, -0.2, 1])
plt.xticks(np.arange(0, 3.5, 0.5))

plt.savefig('Ar_CVV_normalized.png', bbox_inches='tight')
plt.show()

np.savetxt('Ar_CVV_data.txt', average_auto)


# Use time in seconds to fit data
all_steps = np.zeros([len(data)])
for i in range(0, len(average_auto)):
    all_steps[i] = i * DT_fs         # in s

# Fit tail of data to exponential
fit_start = 150
fit_end = len(data)

print(all_steps[fit_start])
print(all_steps[fit_end - 1])
popt_exponential, pcov_exponential = scipy.optimize.curve_fit(exponential, all_steps[fit_start:fit_end],
                                                              average_auto[fit_start:fit_end], p0=[-20000, -3E12])
perr_exponential = np.sqrt(np.diag(pcov_exponential))
print("pre-exponential factor = %0.2f (+/-) %0.2f" % (popt_exponential[0], perr_exponential[0]))
print("rate constant = %0.2f (+/-) %0.2f" % (popt_exponential[1], perr_exponential[1]))

# Plot exponential equation compared to data
plt.figure(figsize=(15, 8))
plt.plot(all_steps[fit_start - 50:fit_end], exponential(all_steps[fit_start - 50:fit_end], popt_exponential[0],
         popt_exponential[1]), 'k--')
plt.plot(all_steps[fit_start - 50:fit_end], average_auto[fit_start - 50:fit_end])
plt.savefig('cvv_fit.png', bbox_inches='tight')

plt.show()


# Extract k and a and convert k to seconds
k = popt_exponential[1]  # 1/s
a = popt_exponential[0]  # m^2


# Integrate from all_steps[150] to infinity
# intergral -> 1/k *exp(kx)

end_step = all_steps[fit_end - 1] + DT_fs
integral_15 = np.sum(average_auto[0:fit_end]) * DT_fs     # m^2/s
integral_15p = 0.0 - (1.0 / k * a * math.exp(k * end_step))   # m^2/ps

print(integral_15 * (100 * 100))
print(integral_15p * (100 * 100))

# Calculate diffusion coefficient
D = (integral_15 + integral_15p) / 3.0  # m^2/s
D = (D) * (100 * 100)                	# cm^2/s
print(D)
