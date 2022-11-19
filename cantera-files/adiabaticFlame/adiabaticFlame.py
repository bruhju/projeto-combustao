import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

gas = ct.Solution('A2highT.yaml')

# T_range = [300, 400, 500, 600, 700, 800, 900, 1000, 1100]
T_range = [814, 707, 1060, 343]

npoints = 150
phi = np.linspace(0, 3, npoints)
tad = np.zeros(npoints)
fuel_species = 'POSF10325'


for T in T_range:
    for i in range(npoints):
        gas.set_equivalence_ratio(phi[i],  'POSF10325',  'O2: 2.0, N2: 7.52')
        gas.TP = T, 101325

        gas.equilibrate('HP',  'auto')

        tad[i] = gas.T-T

    plt.plot(phi, tad)
    # plt.plot(phi, tad, label=(f'{T} K'))


plt.grid(True)
plt.plot(label=r"1 $\phi$")
plt.scatter(x=0.346, y=812.7, c='r', label='812.7 ΔT')
plt.legend(loc='upper right', ncol=2)
# plt.hlines(y=1600, xmin=0, xmax=3, colors='black')
plt.hlines(y=812.7, xmin=0, xmax=0.346, colors='black')

plt.vlines(x=0.346, ymin=0, ymax=812.7, colors='black')


# plt.vlines(x=2.3, ymin=0, ymax=2100, colors='black')

plt.ylim(0, 2000)
plt.xlim(0, 3)

plt.xlabel(r"Equivalence ratio, $\phi$")
plt.ylabel("Temperature [ΔK]")
plt.show()
