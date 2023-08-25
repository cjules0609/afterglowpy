import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

Z = {
    # TopHat, Gaussian, PowerLawCore, GaussianCore, Spherical, PowerLaw
    'jetType':      grb.jet.GaussianCore,
    'specType':     0,
    'thetaObs':     0,                # Viewing angle (rad)
    'E0':           1.5e55,             # Isotropic-equivalent energy (erg)
    'thetaCore':    np.radians(0.8),    # Half-opening angle (rad)
    'thetaWing':    np.radians(np.pi),            # Truncation angle of "wing" (rad)
    'b':            2,                  # Power-law structure index
    'n0':           0.4,                # Number density of ISM (cm-3)
    'p':            2.2,                # Electron distribution power-law index (p>2)
    'epsilon_e':    0.025,               # Thermal energy fraction in electrons
    'epsilon_B':    6e-4,                # Thermal energy fraction in magnetic field
    'xi_N':         1,                # Fraction of electrons accelerated
    'd_L':          2.3e27,             # Luminosity distance (cm)
    'z':            0.151,              # redshift
}

# Space time points geometrically, from 10^3 s to 10^7 s
t = np.geomspace(1, 1.0e4, 1000)

# Calculate flux in a single X-ray band (all times have same frequency)
E = 1e12 # ev
h = 4.136e-15 # ev Hz-1
nu = E/h

# Calculate!
Fnu = grb.fluxDensity(t, nu, **Z)
Fnu_ssc = grb.fluxDensity_ssc(t, nu, **Z)

# Write to a file

print("Writing lc.txt")
with open("lc.txt", 'w') as f:
    f.write("# nu " + str(nu) + '\n')
    f.write("# t(s)     Fnu(mJy)\n")
    for i in range(len(t)):
        f.write("{0:.6e} {1:.6e}\n".format(t[i], Fnu[i]))

print("Plotting")
fig_sync, ax_sync = plt.subplots(1, 1)

ax_sync.plot(t, Fnu)

ax_sync.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$F_\nu$[$10^{26}$ Hz] (mJy)')

fig_sync.tight_layout()
fig_sync.savefig("lc_sync.png")
plt.close(fig_sync)

print("Plotting")
fig_ssc, ax_ssc = plt.subplots(1, 1)

ax_ssc.plot(t, Fnu_ssc)

ax_ssc.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$F_\nu$[$10^{26}$ Hz] (mJy)')

fig_ssc.tight_layout()
fig_ssc.savefig("lc_ssc.png")
plt.close(fig_ssc)

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(t, Fnu, ls='--', label='Synchrotron')
ax.plot(t, Fnu_ssc, ls='dotted', label='SSC')
ax.plot(t, Fnu+Fnu_ssc, ls='-', label='Synchrotron + SSC')

ax.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$F_\nu$[1 TeV] (erg $cm^{-2}$ $s^{-1}$)')

fig.legend()
fig.tight_layout()
fig.savefig("lc_combined.png")
plt.close(fig)
