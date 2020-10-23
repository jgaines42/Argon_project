###################################################
# Plot G(r) using G_r_Ar.txt
############################################

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

# Load g(r) data from simulation
# Column 1: distance in meters
# Column 2: number of points in that bin
data = np.loadtxt('G_r_Ar.txt')

NUM_TIME_STEPS = (50000)/10         # Number of time steps
SIGMA = 0.34E-9             

data[:, 0] = data[:, 0]/SIGMA       # Convert distance to units of sigma 
LENGTH = (3.478E-9)/SIGMA           # Calculate volume in units of sigma
Volume = LENGTH*LENGTH*LENGTH
NATOMS = 864

bins_list = data[:, 0]
results = data[:, 1]
nbins = results.shape[0]
g_of_r = np.zeros(nbins)
scale_factor = np.zeros(500)

# Calculate scale factor
step = bins_list[1]-bins_list[0]
for i in range(1, nbins):
    scale_factor[i] = NATOMS*NATOMS*4*np.pi*step*(i*step)**2

# Scale g(r)
for i in range(1, nbins):
    g_of_r[i] = (results[i])*Volume/(scale_factor[i])/NUM_TIME_STEPS

# Plot
plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(15,8))
plt.plot(bins_list[0:nbins], g_of_r[0:nbins])
plt.plot([0, bins_list[nbins-1]], [1, 1], 'k')
plt.xlabel('$r/\sigma$')
plt.ylabel('$g(r)$')
plt.axis([0, 5, 0, 3])
plt.savefig('Ar_Gr.png', bbox_inches='tight')
g1 = np.concatenate((bins_list, g_of_r))
np.savetxt('Arg_gr_data.txt', g1)
plt.show()
