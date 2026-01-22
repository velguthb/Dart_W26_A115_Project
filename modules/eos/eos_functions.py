import numpy as np


def constants():
    """
    Returns constants 
    [0] : mass of proton 
    [1] : Boltzmann constant 
    [2] : Radiation constant (4 sigmaSB/c_s)
    """ 
    mh = 1.67262192E-24 # mass of proton in g 
    sigma_sb = 	5.670374419E-5 # SB constant: erg⋅cm−2⋅s−1⋅K−4
    kB= 1.380649E-16 # Boltzmann constant: erg/K 
    c_s = 3E10 # speed of light: cm/s 
    a = 4*sigma_sb/c_s 
    return mh, kB, a 

def simple_eos(P, T, mu): 
    """
    Solves the polytropic form of the PV-nKbT equation of state. 
    This function does not handle cases where electron degeneracy 
    or ultra relativistic gas are present. It will best model stars 
    that fall between 0.7Msun <= M 5Msun. 
    :param P: Pressure (g cm^-1 s^-2)
    :param T: Temperature (K)
    :param mu: molar mass (g/moles)
    """
    mh, kB, a = constants()
    pressure = np.asarray(P, dtype=float)
    temp = np.asarray(T, dtype = float)
    P_rad = (1/3)*a*temp**4
    rho = (pressure - P_rad)*(mu*mh)/(kB*temp)
    rho = np.asarray(rho, dtype = float)
    return rho