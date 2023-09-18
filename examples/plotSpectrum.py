import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

Z = {
    # TopHat, Gaussian, PowerLawCore, GaussianCore, Spherical, PowerLaw
    'jetType':       grb.jet.GaussianCore,
    'specType':      0,
    'thetaObs':      0,                # Viewing angle (rad)
    'E0':            1.5e55,             # Isotropic-equivalent energy (erg)
    'thetaCore':     np.radians(0.8),    # Half-opening angle (rad)
    'thetaWing':     np.radians(np.pi),            # Truncation angle of "wing" (rad)
    'b':             2,                  # Power-law structure index
    'n0':            0.4,                # Number density of ISM (cm-3)
    'p':             2.2,                # Electron distribution power-law index (p>2)
    'epsilon_e':     0.025,               # Thermal energy fraction in electrons
    'epsilon_B':     6e-4,                # Thermal energy fraction in magnetic field
    'xi_N':          1,                # Fraction of electrons accelerated
    'd_L':           2.3e27,             # Luminosity distance (cm)
    'z':             0.151,              # redshift

    ######

    'g0':            200,
    # 0: GAMMA_INF, 1: GAMMA_FLAT, 2: GAMMA_EVENMASS, 3: GAMMA_STRUCT
    'gammaType':     1
}


nua = 1.0e7   # Low Frequencies in Hz
nub = 1.0e35  # High Frequencies in Hz

t = 1 * grb.day2sec  # spectrum at 1 day
nu = np.geomspace(nua, nub, num=100)

print("Calculating")
Fnu = grb.fluxDensity(t, nu, **Z)
Fnu_ssc = grb.fluxDensity_ssc(t, nu, **Z)

# mJy Hz to cgs
Fnu = Fnu * 1e-20
Fnu_ssc = Fnu_ssc * 1e-20

# Hz to eV
h = 4.14e-15
ev = nu * h

print("Writing spec.txt")
with open("spec.txt", 'w') as f:
    f.write("# t " + str(t) + ' (s)\n')
    f.write("# nu(Hz)   Fnu(mJy)\n")
    for i in range(len(nu)):
        f.write("{0:.6e} {1:.6e}\n".format(nu[i], Fnu[i]))

print("Plotting synchrotron")
fig_sync, ax_sync = plt.subplots(1, 1)

ax_sync.plot(nu, nu*Fnu)

ax_sync.set_xscale('log')
ax_sync.set_yscale('log')
ax_sync.set_xlabel(r'$\nu$ (Hz)')
ax_sync.set_ylabel(r'$\nu F_\nu$ [1 day] (erg cm$^{-2}$ s$^{-1}$)')

fig_sync.tight_layout()
fig_sync.savefig("spec_sync.pdf")
plt.close(fig_sync)

print("Plotting SSC")
fig_ssc, ax_ssc = plt.subplots(1, 1)

ax_ssc.plot(nu, nu*Fnu_ssc)
# print(nu, Fnu_ssc)

ax_ssc.set_xscale('log')
ax_ssc.set_yscale('log')
ax_ssc.set_xlabel(r'$\nu$ (Hz)')
ax_ssc.set_ylabel(r'$\nu F_\nu$ [1 day] (erg cm$^{-2}$ s$^{-1}$)')

fig_ssc.tight_layout()
fig_ssc.savefig("spec_ssc.pdf")
plt.close(fig_ssc)

print("Plotting combined")
fig, ax = plt.subplots(1, 1)

ax.plot(ev, nu*Fnu, ls='--', label='Synchrotron')
ax.plot(ev, nu*Fnu_ssc, ls='-', label='SSC')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'eV')
ax.set_ylabel(r'$\nu F_\nu$ [1 day] (erg cm$^{-2}$ s$^{-1}$)')
ax.set_xlim([1e-7,1e20])
ax.set_ylim([1e-10, 2e-4])

fig.legend(borderaxespad=1.5)
fig.tight_layout()
# plt.show()
fig.savefig("spec_combined.pdf")
plt.close(fig)
