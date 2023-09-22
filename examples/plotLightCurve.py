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

    'g0':            500,
    # 0: GAMMA_INF, 1: GAMMA_FLAT, 2: GAMMA_EVENMASS, 3: GAMMA_STRUCT
    'gammaType':     1
}

# Space time points geometrically, from 10^3 s to 10^7 s
t = np.geomspace(1, 1.0e8, 500)

# Calculate flux in a single X-ray band (all times have same frequency)
E = 1e3 # ev
h = 4.136e-15 # ev Hz-1
nu = E/h

# Calculate!
Z['radType'] = 0
Fnu = grb.fluxDensity(t, nu, **Z)
Z['radType'] = 1
Fnu_ssc = grb.fluxDensity(t, nu, **Z)

Y = Fnu_ssc / Fnu
Fnu_KN = Fnu * (1 + Y)

# mJy Hz to cgs
nuFnu = Fnu * nu * 1e-20
nuFnu_ssc = Fnu_ssc * nu * 1e-20
nuFnu_KN = Fnu_KN * nu * 1e-20

# Write to a file

print("Writing lc.txt")
with open("lc.txt", 'w') as f:
    f.write("# nu " + str(nu) + '\n')
    f.write("# t(s)     Fnu(mJy)\n")
    for i in range(len(t)):
        f.write("{0:.6e} {1:.6e}\n".format(t[i], Fnu[i]))

print("Plotting")
fig_sync, ax_sync = plt.subplots(1, 1)

ax_sync.plot(t, nuFnu)

ax_sync.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$\nu F_\nu$ [1 TeV] (erg cm$^{-2}$ s$^{-1}$)')
ax_sync.set_xlim([1,1e8])

fig_sync.tight_layout()
fig_sync.savefig("lc_sync.pdf")
plt.close(fig_sync)

print("Plotting")
fig_ssc, ax_ssc = plt.subplots(1, 1)

ax_ssc.plot(t, nuFnu_ssc)

ax_ssc.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$\nu F_\nu$ [1 TeV] (erg cm$^{-2}$ s$^{-1}$)')
ax_ssc.set_xlim([1,1e8])

fig_ssc.tight_layout()
fig_ssc.savefig("lc_ssc.pdf")
plt.close(fig_ssc)

print("Plotting")
fig, ax = plt.subplots(1, 1)

ax.plot(t, nuFnu, ls='--', label='Synchrotron')
ax.plot(t, nuFnu_ssc, ls='dotted', label='SSC')
ax.plot(t, nuFnu+nuFnu_ssc, ls='-', label='Synchrotron + SSC (w/o KN)')
ax.plot(t, nuFnu_KN, ls = '-.', label='Synchrotron + SSC (w/ KN)')

ax.set(xscale='log', xlabel=r'$t$ (s)',
       yscale='log', ylabel=r'$\nu F_\nu$ [1 TeV] (erg cm$^{-2}$ s$^{-1}$)')
ax.set_xlim([1,1e8])

fig.legend(loc='lower left', bbox_to_anchor=[0.15,0.15])
fig.tight_layout()
fig.savefig("lc_combined.pdf")
plt.close(fig)
