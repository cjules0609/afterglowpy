import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb

# lc(): calculates synchrotron and ssc light curves
# Inputs: 
    # Z (model parameters)
    # var (str: variable or nu to iterate over)
    # val (list: values for variable to iterate over)
    # t (list: time grid)
    # nu (tuple: frequency bounds)
# Returns: 
    # F (flux [cgs])
def lc(Z, var, val, t=np.logspace(2, 9, num=500), nu=[1e17, 1e18]):
    print('Calculating light curve for ' + var + ' = ' + str(round(val, 2)))
    if var == 'nu':
        F = grb.fluxDensity(t, val, **Z)
        F, nu_cgs = convert_cgs(F, val)
        F *= nu_cgs

    else:
        n = 5
        nu_list = np.linspace(nu[0], nu[1], n)

        F_list = []
        Z[var] = val
        # integrate over nu range
        for j, nu in enumerate(nu_list):
            F = grb.fluxDensity(t, nu, **Z)
            F, nu_cgs = convert_cgs(F, nu)
            F_list.append(F * nu_cgs)
        F = np.mean(F_list, axis=0)

    return F

# spec(): calculates synchrotron and ssc spectra
# Inputs: 
    # Z (model parameters)
    # var (str: variable or nu to iterate over)
    # val (list: values for variable to iterate over)
    # t (float: time)
    # nu (list: frequency grid)
# Returns: 
    # F (flux [cgs])
def spec(Z, var, val, nu=np.geomspace(1e9, 1e20, num=500), t=1e2):
    print('Calculating spectrum for ' + var + ' = ' + str("_{:.2e}".format(val)))
    ssc = True
    if var == 't':
        Z['radType'] = 0
        F = grb.fluxDensity(val, nu, **Z)
        F, nu_cgs = convert_cgs(F, val)
        if ssc:
            Z['radType'] = 1
            F_ssc = grb.fluxDensity(val, nu, **Z)
            F_ssc, nu_cgs = convert_cgs(F_ssc, nu)      
            F = [sum(x) for x in zip(F, F_ssc)]

    else:
        Z[var] = val
        F = grb.fluxDensity(t, nu, **Z)
        F, nu_cgs = convert_cgs(F, nu)
        if ssc:
            Z['radType'] = 1
            F_ssc = grb.fluxDensity(t, nu, **Z)
            F_ssc, nu_cgs = convert_cgs(F, nu)       
            F = [sum(x) for x in zip(F, F_ssc)]

    return F * nu

# convert_cgs(): converts afterglowpy outputs to cgs
# Inputs: 
    # F (list: flux [mJy])
# Returns: 
    # F (list: flux [erg/s/cm^2/Hz])
def convert_cgs(F, nu):
    # mJy to erg/s/cm^2/Hz
    F *= 1e-26

    # Hz to eV
    h = 4.14e-15
    ev = nu * h
    # Stay in nu (or check if t and convert?)
    x = nu
    return F, x