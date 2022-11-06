import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

ct.suppress_thermo_warnings()
T_range = [300, 400, 500, 600, 700, 800, 900, 1000, 1100]
P = ct.one_atm  # type: ignore

# phases
gas = ct.Solution('jetA1.yaml')   # type: ignore
carbon = ct.Solution('graphite.yaml')  # type: ignore

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0), (carbon, 0.0)]


# gaseous fuel species
fuel_species = 'CH4'

# equivalence ratio range
npoints = 50
phi = np.linspace(0, 3.5, npoints)

##############################################################################

mix = ct.Mixture(mix_phases)  # type: ignore

# create some arrays to hold the data
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species, npoints))

for T in T_range:
    for i in range(npoints):
        # set the gas state
        gas.set_equivalence_ratio(phi[i], fuel_species, 'O2:1.0, N2:3.76')

        # create a mixture of 1 mole of gas, and 0 moles of solid carbon.
        mix = ct.Mixture(mix_phases)  # type: ignore
        mix.T = T
        mix.P = P

        # equilibrate the mixture adiabatically at constant P
        mix.equilibrate('HP', max_steps=1000)
        tad[i] = mix.T
        xeq[:, i] = mix.species_moles

    plt.plot(phi, tad, label=(f'{T} K'))


plt.grid(True)
plt.legend(loc='upper right', ncol=3)
plt.xlabel(r"Equivalence ratio, $\phi$")
plt.xlim(0, 2.8)
plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 2.8])
plt.ylabel("Temperature [K]")
plt.ylim(300, 2700)
plt.yticks([300, 400.,  500.,  600.,  700.,  800.,  900., 1000., 1100., 1200.,
            1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2100.,
            2200., 2300., 2400., 2500., 2600., 2700., 2800.])
plt.show()
