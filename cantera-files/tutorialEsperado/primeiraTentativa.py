import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.dpi"] = 120
ct.suppress_thermo_warnings()


print('Start procediment...\n\n')

# Get all of the Species objects defined in the GRI 3.0 mechanism
species = {S.name: S for S in ct.Species.list_from_file("jetA1.yaml")}
air = "O2:0.21,N2:0.79"
jet_fuel_composition = "nC3H7:0.268,iC3H7:0.397,c-C4H5:0.201,aC3H4:0.134"

# Create an IdealGas object with species representing complete combustion
complete_species = [species[S] for S in species]
jet_fuel = ct.Solution(thermo="IdealGas", species=complete_species)


phi = np.linspace(0.5, 2.0, 100)
T_complete1 = np.zeros(phi.shape)
for i in range(len(phi)):
    jet_fuel.TP = 300, ct.one_atm
    jet_fuel.set_equivalence_ratio(phi[i], "POSF10264", air)
    jet_fuel.equilibrate("HP")
    T_complete1[i] = jet_fuel.T


plt.plot(phi, T_complete1, label="complete combustion", lw=2)
plt.grid(True)
plt.xlabel(r"Equivalence ratio, $\phi$")
plt.ylabel("Temperature [K]")

plt.show()
