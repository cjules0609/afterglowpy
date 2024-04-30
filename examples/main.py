import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb
import pandas as pd
import csv

import calc_flux

h = 4.136e-15 # ev Hz-1
cwd = '/Volumes/T7Shield/ncfa/afterglowpy/examples/'

Z = {'jetType':     grb.jet.Gaussian,     # Top-Hat jet
     'specType':    0,                  # Basic Synchrotron Emission Spectrum

     'thetaObs':    0.05,   # Viewing angle in radians
     'E0':          1.0e53, # Isotropic-equivalent energy in erg
     'thetaCore':   0.1,    # Half-opening angle in radians
     'thetaWing':   0.5,
     'n0':          1.0,    # circumburst density in cm^{-3}
     'p':           2.2,    # electron energy distribution index
     'epsilon_e':   0.1,    # epsilon_e
     'epsilon_B':   0.01,   # epsilon_B
     'xi_N':        1.0,    # Fraction of electrons accelerated
     'd_L':         1.0e28, # Luminosity distance in cm
     'z':           0.55}   # redshift

model_170817 = {
        # TopHat, Gaussian, PowerLawCore, GaussianCore, Spherical, PowerLaw
        'jetType':       grb.jet.Gaussian,
        'specType':      0,
        'thetaObs':      0,                # Viewing angle (rad)
        'E0':            1e53,             # Isotropic-equivalent energy (erg)
        'thetaCore':     0.088,    # Half-opening angle (rad)
        'thetaWing':     1,            # Truncation angle of "wing" (rad)
        'n0':            1*10**-1.7,                # Number density of ISM (cm-3)
        'p':             2.139,                # Electron distribution power-law index (p>2)
        'epsilon_e':     1e-2,               # Thermal energy fraction in electrons
        'epsilon_B':     1*10**-3.7,                # Thermal energy fraction in magnetic field
        'xi_N':          1.0,                # Fraction of electrons accelerated
        'd_L':           1.23e26,             # Luminosity distance (cm)
        'z':             0.01,

        # 0: GAMMA_INF, 1: GAMMA_FLAT, 2: GAMMA_EVENMASS, 3: GAMMA_STRUCT
        'gammaType':     0,
        # 'g0':            3000,
        'spread':        False,
        # 0: Synchrotron 1: SSC
        'radType':       0,
        # 0: no cooling 1: SSC cooling 2: SSC cooling w/ KN
        # 'coolType':   2
    }

model_240422 = {
        'jetType':       grb.jet.Gaussian,
        'specType':      0,
        'thetaObs':      0.54,                # Viewing angle (rad)
        'E0':            1e51,             # Isotropic-equivalent energy (erg)
        'thetaCore':     0.88,    # Half-opening angle (rad)
        'thetaWing':     1,            # Truncation angle of "wing" (rad)
        'n0':            1,                # Number density of ISM (cm-3)
        'p':             2.139,                # Electron distribution power-law index (p>2)
        'epsilon_e':     1e-1,               # Thermal energy fraction in electrons
        'epsilon_B':     1e-2,                # Thermal energy fraction in magnetic field
        'xi_N':          1.0,                # Fraction of electrons accelerated
        'd_L':           6.6e26,             # Luminosity distance (cm)
        'z':             0.01,

        # 0: GAMMA_INF, 1: GAMMA_FLAT, 2: GAMMA_EVENMASS, 3: GAMMA_STRUCT
        'gammaType':     0,
        # 'g0':            3000,
        'spread':        False,
        # 0: Synchrotron 1: SSC
        'radType':       0,
        # 0: no cooling 1: SSC cooling 2: SSC cooling w/ KN
        # 'coolType':   2
    }

Z = model_240422

# main(): initializes the program
# Inputs: 
    # None
# Returns: 
    # None
def main():
    # for user-friendly UI
    # init_ui()

    # for backend use and debugging
    spec_lc = dev_mode()

    # plot observed data
    plot_obs(spec_lc)

    format_plot(spec_lc)

# dev_mode(): enables backend features
# Inputs: 
    # None
# Returns: 
    # None
def dev_mode():
    print('-------------afterglowpy (dev mode)-------------')
    spec_lc = 'lc'

    # input 'parameter' to vary
    # or 'nu' for lc
    # or 't' for spec
    var = 'thetaObs'

    # values to iterate
    # vals = [1e2, 1e4, 1e6, 1e8]
    vals = [np.radians(x) for x in [0, 20, 40, 60, 80]]
    # vals = [np.radians(x) for x in [1, 5, 10, 20, 30, 40]]

    file_list = model_to_datafile(Z, spec_lc, var=var, vals=vals)
    # file_list = ['lc_thetaCore_0.028_E54', 'lc_thetaCore_0.088_E53', 'lc_thetaCore_0.278_E52', 'lc_thetaCore_0.54', 'lc_thetaCore_0.88_E51']
    fig, ax = datafile_to_plot(file_list)
    
    return spec_lc

# init_ui(): Initializes UI for user-friendly operation
# Inputs: 
    # None
# Returns: 
    # None
def init_ui():
    print('-------------afterglowpy UI-------------')

    while True:
        operation = input("What do you want to calculate?\n0: spectra\n1: lightcurves\n-> ")
        if operation == '0':
            spec_lc = 'spec'
            var = 't'
            vals = np.geomspace(1, 1e6, num=5)
            break
        if operation == '1':
            spec_lc = 'lc'
            var = 'nu'
            vals = np.logspace(6, 28, num=5)
            break
        print('Please input 0 or 1')
        
    
    file_list = model_to_datafile(Z, spec_lc, var=var, vals=vals)
    while True:
        plot_data = input('Plot data?\n0: no\n1: yes\n-> ')
        if plot_data == '0' or plot_data == '1':
            break
        print('Please input 0 or 1')

    if plot_data:
        datafile_to_plot(file_list)

    print('All done. Goodbye...')
    return

# model_to_default(): calculates flux light curves and spectra and saves data files
# Inputs: 
    # Z (model parameters)
    # spec_lc ('spec' for spectrum or 'lc' for lc)
# Returns: 
    # file_list (list of data files)
def model_to_datafile(Z, spec_lc, var=None, vals=None):
    path = cwd + 'data/'
    file_list = []
    # lc settings
    if spec_lc == 'lc':
        for i, val in enumerate(vals):
            t = np.logspace(2, 9, num=500)
            nua = 0.5e3 / h
            nub = 10e3 / h
            nu = [nua, nub]
            F = calc_flux.lc(Z, var, val, t=t, nu=nu)

            # convert radians to degrees for filename
            if var == 'thetaObs' or var == 'thetaCore' or var == 'thetaWing':
                val = val * 360 / 2 / np.pi
            filename =  'lc_' + var + str("_{:.1e}".format(val))
            print('Saving data...')
            write_csv(t, F, path, filename, xname='t')
            file_list.append(filename)

    # spectrum settings
    elif spec_lc == 'spec':
        for i, val in enumerate(vals):
            nu = np.geomspace(1e6, 1e30, num=1000)
            F = calc_flux.spec(Z, var, val, nu=nu, t=10)

            # convert radians to degrees for filename
            if var == 'thetaObs' or var == 'thetaCore' or var == 'thetaWing':
                val = val * 360 / 2 / np.pi
            filename =  'spec_' + var + str("_{:.1e}".format(val))
            print('Saving data...')
            write_csv(nu, F, path, filename, xname='nu')
            file_list.append(filename)
    
    else:
        return None

    return file_list

# datafile_to_plot(): reads data files and plots data
# Inputs: 
    # file_list (list of filenames)
# Returns: 
    # fig, ax of plot
def datafile_to_plot(file_list):
    print('Plotting data...')

    key = {
        'thetaObs':     r'$\theta_v$',
        'E0':           r'$E_0$',
        'thetaCore':    r'$\theta_c$',
        'thetaWing':    r'$\theta_w$',
        'n0':           r'$n_0$',
        'p':            r'$p',
        'epsilon_e':    r'$\epsilon_e$',
        'epsilon_B':    r'$\epsilon_B$',
        'xi_N':         r'$Xi_n$',
        'd_L':          r'$\d_L$',
        'z':            r'$z$',
    }

    # generate plot
    fig, ax = plt.subplots()
    colors = plt.cm.plasma(np.linspace(0,1,len(file_list)))

    # plot model data
    path = cwd + 'data/'
    for i, file in enumerate(file_list):
        label = key[file.split('_')[1]] + '=' + str(float(file.split('_')[2])) if file.split('_')[1] in key else file
        model_data = import_data(path, file)
        if model_data.columns[0] == 't':
            model_data['t'] = model_data['t']/86400
        model_data.plot(ax=ax, x=model_data.columns[0], y='flux', logx=True, logy=True, label=label, color=colors[i])

    return fig, ax

# plot_obs(): plot any observed data
# Inputs: 
    # None
# Returns: 
    # None
def plot_obs(spec_lc):
    spec_lc = None

    ax = plt.gca()
    # plot 170817 data
    # lc
    # path = '/Volumes/T7Shield/ncfa/GRB170817A/'
    # files = ['Piro+19', 'Troja+18', 'Troja+19']
    # obs_data = import_data(path, files)
    # obs_data.plot.scatter(ax=ax, x='t', y='flux', yerr=[obs_data['err-'], obs_data['err+']], s=10, label='170817')
    
    # plot EP data
    t_ep = [2.8e5/86400, 2.9e5/86400]
    nu_ep = [1e3/h, 1e3/h]

    # plot yahpt data
    t_yahpt = [2.272e5/86400, 2.275e5/86400]
    nu_yahpt = [4.56e14, 3.72e14]
    f_yahpt = [4.28e-14, 3.43e-14]

    if spec_lc == 'spec':
        ax.scatter(nu_ep, [4e-13, 3e-13], label='EP', s=10, c='k')
        ax.scatter(nu_yahpt, f_yahpt, label='YAHPT', s=10, c='r', marker='v')
    
    if spec_lc == 'lc':
        # ax.scatter(t_ep, [4e-13, 3e-13], label='EP', s=10, c='k')
        ax.scatter(t_yahpt, f_yahpt, label='YAHPT', s=10, c='r', marker='v')

    # ax.fill_between([1, 1e4], [1e-13, 1e-13], [1, 1], label='EP', facecolor="none", hatch="X", edgecolor="k", linewidth=0.5)
    ax.scatter(1, 1e-13, label='EP', marker='v', c='k')


# import_data(): reads data files
# Inputs: 
    # path (directory of files)
    # files (list of filenames)
# Returns: 
    # None
def import_data(path, files):
    df = pd.DataFrame()
    if type(files) != list:
        files = [files]
    for filename in files:
        filedata = pd.read_csv(path + filename + '.csv')
        df = pd.concat([df, filedata])

    return df

# write_csv(): writes data to .csv files
# Inputs: 
    # x (list: x variable | t for lc | nu for spec)
    # F (list: flux)
    # path (str: file directory)
    # filename (str: filename to save)
# Returns: 
    # None
def write_csv(x, F, path, filename, xname=None):
    with open(path + filename+'.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow([xname, 'flux'])
        for i, val in enumerate(F):
            writer.writerow([x[i], F[i]])

def format_plot(plotname):
    # format plot
    if plotname == 'lc':
        # lc
        plt.xlabel(r'$T - T_0$ [d]'); plt.ylabel(r'Flux ($5-10$ keV) [erg/s/cm$^2$]')
        plt.xlim([1e-2,1e4]); plt.ylim([1e-22, 1e-10])
    elif plotname == 'spec':
        # spec
        plt.xlabel(r'$\nu$ [Hz]'); plt.ylabel(r'$F$ [erg/s/cm$^2$]')
        plt.xlim([1e7,1e30]); plt.ylim([1e-20, 1e-12])

    plt.legend()
    plt.tight_layout()
    plt.savefig(cwd + 'out/' + plotname + '.pdf')
    plt.show()

if __name__ == '__main__':
    main()