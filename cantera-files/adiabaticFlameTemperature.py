import cantera as ct
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T_range = [300,400,500,600,700,800,900,1000]
P = ct.one_atm

# phases
gas = ct.Solution('jetA1.yaml')
carbon = ct.Solution('graphite.yaml')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0)]

# gaseous fuel species
fuel_species = 'POSF10264'

# equivalence ratio range
npoints = 50
phi = np.linspace(0, 3.5, npoints)

##############################################################################

mix = ct.Mixture(mix_phases)

# create some arrays to hold the data
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species, npoints))

for T in T_range:
    for i in range(npoints):
        # set the gas state
        gas.set_equivalence_ratio(phi[i], fuel_species, 'O2:1.0, N2:3.76')

        # create a mixture of 1 mole of gas, and 0 moles of solid carbon.
        mix = ct.Mixture(mix_phases)
        mix.T = T
        mix.P = P

        # equilibrate the mixture adiabatically at constant P
        mix.equilibrate('HP', solver='vcs')
        tad[i] = mix.T
        xeq[:, i] = mix.species_moles
        
    plt.plot(phi, tad,label=T)


plt.grid(True)
plt.legend(loc='upper right',ncol=2)
plt.xlabel(r"Equivalence ratio, $\phi$")
plt.ylabel("Temperature [K]")
plt.show()



    
