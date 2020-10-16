import numpy as np
import math
import matplotlib.pyplot as plt
import time

############################
# set_up_coordinates()
# Input:
#   L: Length of box
#   N: number of atoms
#   sigma: radii * 2
# returns: coordinates
#############################
def set_up_coordinates(L,N, sigma):
    FCC_lattice = np.zeros([N,3])
    r = sigma/math.sqrt(2)*1.15
    vertical = L/6 # We know that we need to fit 6 replicates in each dimension

    # repeating pattern of atoms
    FCC_atoms = np.asarray([(0, 0, 0), (r, r, 0), (r, 0, r), (0, r, r)])
    n = 0;
    for x in range(0,6):
        for y in range(0, 6):
            for z in range(0, 6):
                for i in range(0,4):
                    coordinatestranslation = vertical*np.asarray([x, y, z]);
                    FCC_lattice[n,:] = coordinatestranslation + FCC_atoms[i,:]
                    n = n+1
    np.savetxt('inital_coordinates.out', FCC_lattice, delimiter=',')

    return FCC_lattice;

def random_three_vector():
    """
    https://gist.github.com/andrewbolster/10274979
    Generates a random 3D unit vector (direction) with a uniform spherical distribution
    Algo from http://stackoverflow.com/questions/5408276/python-uniform-spherical-distribution
    :return:
    """
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return (x,y,z)


############################
# set_up_velocities()
# Input:
#   L
#   N
# returns: velocities
#############################
def set_up_velocities(N, velMag):

    # keep track of global velocity to make sure we're not shifting
    vSum = np.zeros([3])
    vSum2 = np.zeros([3])
    initial_velocities = np.zeros([N,3])

    # loop over all atoms and give them a random velocity
    for i in range(0, N):
        rand_direction = random_three_vector()
        initial_velocities[i,:] = np.multiply(rand_direction, velMag)
        vSum = vSum + initial_velocities[i,:]
    print(vSum)
    # Make it so we're not shifting
    for i in range(0, N):
        initial_velocities[i,:]= initial_velocities[i,:] + np.multiply((-1.0/ N),  vSum)
        vSum2 = vSum2 +initial_velocities[i,:]
    #print(vSum2)

    np.savetxt('inital_velocities.out', initial_velocities, delimiter=',',fmt='%2.3e')
    return initial_velocities

#############################
# calculate_forces()
# Input:
#   current_coordinates
#   sigma
#   epsilon
#   L
# Output:
#   forces (Nx3)
##############################
def calculate_forces(N, sigma, current_coordinates, epsilon,L):

    #force_cutoff = 100.0*sigma#2.25*sigma

    V_tot = 0
    new_forces = np.zeros([N,3])
    # loop over all pairs of atoms
    for i in range(0,N-1): #N-1
        for j in range(i+1,N): # i+1,N
            rij = current_coordinates[i,:]-current_coordinates[j,:]

            # Deal with periodic
            if rij[0] > L/2.0:
                rij[0] = -(L-rij[0])
            elif rij[0] < -L/2.0:
                rij[0] = (L+rij[0])
            if rij[1] > L/2.0:
                rij[1] = -(L-rij[1])
            elif rij[1] < -L/2.0:
                rij[1] = (L+rij[1])
            if rij[2] > L/2.0:
                rij[2] = -(L-rij[2])
            elif rij[2] < -L/2.0:
                rij[2] = (L+rij[2])

            dist_rij = math.sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])



            if dist_rij < 10000000:#force_cutoff:
                s_div_rij = sigma/dist_rij;
                sr_6 = s_div_rij**6
                sr_12 = s_div_rij**12
                VIJ=4*epsilon*(sr_12-sr_6)
                V_tot=V_tot+VIJ

                this_force=-24*epsilon*(-2.0*sr_12+sr_6)/dist_rij

                unit_rij = np.divide(rij,dist_rij);
                new_forces[i,:]= new_forces[i,:] + np.multiply(this_force,unit_rij);
                new_forces[j,:] = new_forces[j,:] - np.multiply(this_force,unit_rij);

    return new_forces, V_tot

def calculate_forces_new(N, sigma, current_coordinates, epsilon,L):

    force_cutoff =2.5*sigma
    force_cutoff_2 = force_cutoff*force_cutoff
    V_tot = 0
    new_forces = np.zeros([N,3])
    # loop over all pairs of atoms
    for i in range(0,N-1): #N-1
        all_rij = current_coordinates[i,:]-current_coordinates[i+1:N,:]
        ind0 = all_rij> L/2.0
        all_rij[ind0] = -(L-all_rij[ind0])
        ind1 = all_rij < -L/2.0
        all_rij[ind1] = (L+all_rij[ind1])

        for j in range(0,N-i-1): # i+1,N
            # rij = current_coordinates[i,:]-current_coordinates[j,:]
            #
            # # Deal with periodic
            # if rij[0] > L/2.0:
            #     rij[0] = -(L-rij[0])
            # elif rij[0] < -L/2.0:
            #     rij[0] = (L+rij[0])
            # if rij[1] > L/2.0:
            #     rij[1] = -(L-rij[1])
            # elif rij[1] < -L/2.0:
            #     rij[1] = (L+rij[1])
            # if rij[2] > L/2.0:
            #     rij[2] = -(L-rij[2])
            # elif rij[2] < -L/2.0:
            #     rij[2] = (L+rij[2])
            rij = all_rij[j,:]
            dist_rij_2 = (rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])



            if dist_rij_2 < 10000000:#force_cutoff_2:
                dist_rij = math.sqrt(dist_rij_2)
                s_div_rij = sigma/dist_rij;
                sr_6 = s_div_rij**6
                sr_12 = s_div_rij**12
                VIJ=4*epsilon*(sr_12-sr_6)
                V_tot=V_tot+VIJ

                this_force=-24*epsilon*(-2.0*sr_12+sr_6)/dist_rij

                unit_rij = np.divide(rij,dist_rij);
                new_forces[i,:]= new_forces[i,:] + np.multiply(this_force,unit_rij);
                new_forces[j+i+1,:] = new_forces[j+i+1,:] - np.multiply(this_force,unit_rij);

    return new_forces, V_tot



def calculate_KE(N, current_velocities, massAr):
    # squared_velocities = np.square(current_velocities)
    # sumvsq = np.sum(squared_velocities)
    sumvsq = 0;
    for i in range(0,N):
        sumvsq=sumvsq+current_velocities[i,0]**2 + current_velocities[i,1]**2 + current_velocities[i,2]**2
    KE=0.5*massAr*sumvsq
    return KE

def calculate_temp(KE, Boltz):
    temp=(KE*2.0)/(3.0*Boltz*N)
    return temp

##############################
# update_coordinates()
# Input:
#   current_coordinates
#   current_velocities
#   forces
#   L
# Output:
#   new_coordinates
#   new_velocities
##############################
def update_velocity( current_accelerations, change_acceleration, current_velocities, N,  dt, is_equi, Tref, time, all_temp):

    # update accelerations
    current_accelerations= change_acceleration# + current_accelerations
    current_velocities = np.add(current_velocities , np.multiply(current_accelerations,dt))

    return current_velocities

def scale_velocity(current_velocities, N,  dt, is_equi, Tref, time, all_temp, Boltz):

    KE = calculate_KE(N, current_velocities, massAr)
    temp=calculate_temp(KE, Boltz)
    all_temp[time,0] = temp
    if temp > (Tref+2.5) or temp < (Tref-2.5):
        current_velocities = np.multiply(current_velocities, math.sqrt(Tref/temp))

    return current_velocities, KE


def update_coordinates(current_coordinates, current_velocities, L, dt):
    change_coordinates = current_velocities*dt
    current_coordinates = current_coordinates+change_coordinates;

    ind0 = current_coordinates > L
    current_coordinates[ind0] = current_coordinates[ind0]-L
    ind1 = current_coordinates < 0
    current_coordinates[ind1] = current_coordinates[ind1]+L
    #current_coordinates = np.add(current_coordinates , np.multiply(current_velocities,dt))
    # for i in range(0,N):
    #     #Deal with periodic boundries
    #     for j in range(0,3):
    #         if current_coordinates[i,j]> L :
    #             current_coordinates[i,j] = current_coordinates[i,j]-L
    #         if current_coordinates[i,j] < 0:
    #             current_coordinates[i,j] = current_coordinates[i,j] + L
    return current_coordinates


## main function
tic = time.perf_counter()
# Get parameters
NA=6.02214076E23                    # Avogadros number
L = 3.478E-9                         # Length in m
N = 864
Boltz=1.38064852E-23                # in J/K

epsilon = 120.0*Boltz                 # in J
sigma = 3.4E-10;                    # sum of radii in m
#nstep=1000000

massAr=0.039948/NA                  # in kg
Tref=94.400000                      # in Kelvin
NDIM = 3;
temperature = Tref

velMag = math.sqrt(3.0*Boltz*Tref/massAr)                 # in m/s : ((3*boltzmann_constant*Temperature)/ mass)^0.5
dt = 10.0E-15;                      # in seconds
num_runs = 100

# seed random nubmer
np.random.seed(1)

# set up inital  positions
current_coordinates = set_up_coordinates(L,N, sigma)
current_velocities = set_up_velocities(N, velMag)
current_accelerations = np.zeros([N,3])

KE = calculate_KE(N, current_velocities, massAr)
temp=calculate_temp(KE, Boltz)

## Save initial coordinates
pdb = open('Initial_coordinates.pdb', 'w')
pdb.write('%3i \n'%N)
for i in range(0,N):
    pdb.write('ATOM   %4i  Ar  ArA A %3i      '%(i,i))
    pdb.write('%6.3f  %6.3f  %6.3f \n'%(current_coordinates[i,0]*1E10, current_coordinates[i,1]*1E10,current_coordinates[i,2]*1E10 ))
pdb.close()

## Save initial velocities
pdb = open('Initial_velocities.xyz', 'w')
pdb.write('%3i \n'%N)
pdb.write('Comment:Argon 864 cluster \n')
for i in range(0,N):
    #pdb.write()
    #pdb.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f} \n".format('ATOM', i, 'Ar', 'A' 'ArA', 'A', i, ' ',current_velocities[i,0], current_velocities[i,1],current_velocities[i,2] ))
    #pdb.write('ATOM   %4i  Ar  ArA A %3i     '%(i,i))
    pdb.write('Ar %07.3f %07.3f %07.3f \n'%(current_velocities[i,0], current_velocities[i,1],current_velocities[i,2]))
pdb.close()

# open output file
f = open('Temp.txt', 'w')
all_energy = np.zeros([num_runs,3])


forces, PE = calculate_forces(N, sigma, current_coordinates, epsilon, L)
all_energy[0,0] = KE
all_energy[0,1] = PE
all_energy[0,2] = PE + KE
np.savetxt('initial_forces.out', forces, delimiter=',',fmt='%2.3e')
np.savetxt('initial_velocity.out', current_velocities, delimiter=',',fmt='%2.3e')
np.savetxt('initial_coordinates.out', current_coordinates, delimiter=',',fmt='%2.3e')
np.savetxt('initial_accelerations.out', current_accelerations, delimiter=',',fmt='%2.3e')

is_equi = 0
all_temp = np.zeros([num_runs,1])
all_temp[0] = temp
## Loop over time
for t in range(1,num_runs):
    # calculate forces
    forces, PE = calculate_forces_new(N, sigma, current_coordinates, epsilon, L)
    change_acceleration = np.divide(forces,massAr);
    # update coordinates
    current_velocities = update_velocity(current_accelerations, change_acceleration, current_velocities, N,  dt, is_equi, Tref, t, all_temp)
    KE = calculate_KE(N, current_velocities, massAr)
    temp=calculate_temp(KE, Boltz)
    all_temp[t,0] = temp

    #current_velocities, KE = scale_velocity(current_velocities, N,  dt, is_equi, Tref, t, all_temp, Boltz)

    current_coordinates = update_coordinates(current_coordinates, current_velocities, L, dt)
    all_energy[t,0] = KE
    all_energy[t,1] = PE
    all_energy[t,2] = PE + KE
toc = time.perf_counter()
print(f"Ran in {toc - tic:0.4f} seconds")

np.savetxt('final_forces.out', forces, delimiter=',',fmt='%2.3e')
np.savetxt('final_velocity.out', current_velocities, delimiter=',',fmt='%2.3e')
np.savetxt('final_coordinates.out', current_coordinates, delimiter=',',fmt='%2.3e')
np.savetxt('final_accelerations.out', current_accelerations, delimiter=',',fmt='%2.3e')
np.savetxt('All_energy.txt', all_energy, fmt='%10.8e')

f.close()
np.savetxt('all_temp.out', all_temp, delimiter=',',fmt='%4.2f')

plt.plot(range(1,num_runs+1), all_temp[0:num_runs,0])
plt.xlabel('time')
plt.ylabel('temp')
plt.show()

plt.plot(range(1,num_runs+1), all_energy[0:num_runs,0])
plt.xlabel('time steps (1E-14)')
plt.ylabel('Kinetic energy (J)')
plt.show()

plt.plot(range(1,num_runs+1), all_energy[0:num_runs,1])

plt.xlabel('time')
plt.ylabel('Potential energy (J)')
plt.show()

plt.plot(range(1,num_runs+1), all_energy[0:num_runs,2])
plt.plot(range(1,num_runs+1), all_energy[0:num_runs,1])
plt.plot(range(1,num_runs+1), all_energy[0:num_runs,0])
plt.xlabel('time steps (1E-14)')
plt.ylabel('Total energy (J)')
plt.show()
