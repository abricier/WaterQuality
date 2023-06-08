# This module contains the functions that will be used to invert from reflectances to water quality parameters
# The dicionary `functions`` will be used by the waterquality module:
# function = {
#   'name_of_the_parameter': {'function': any_function}, 'units': 'units to be displayed in the report',
#   'name_of_the_parameter2': .....
# }

# Any bands can be used to compute the final value. The name of the band must match the internal name used by WaterDetect
# It is enough to put the band name as an argument in the function
# Available bands in Sentinel2 are:
# Blue, Green, Red, RedEdg1, RedEdg2, RedEdg3, Nir, Nir2, Mir, Mir2

import numpy as np

# Below is an example extracted from Nechad et al. (2010)
def nechad(Red, a=610.94, c=0.2324):
    spm = a * Red / (1 - (Red / c))
    return spm

# SEN3R SPM
def _spm_modis(Nir, Red):
    return 759.12 * ((Nir / Red) ** 1.92)

def _power(x, a, b, c):
    return a * (x) ** (b) + c

# Gitelson
def chl_gitelson(Red, RedEdg1, RedEdg2):
    chl = 23.1 + 117.4 * (1 / Red - 1 / RedEdg1) * RedEdg2
    return chl


# Gitelson and Kondratyev
def chl_gitelson2(Red, RedEdg1):
    chl = 61.324 * (RedEdg1 / Red) - 37.94
    return chl


# Turbidity (FNU) Dogliotti
def turb_dogliotti(Red, Nir2):
    """Switching semi-analytical-algorithm computes turbidity from red and NIR band

    following Dogliotti et al., 2015
    :param water_mask: mask with the water pixels (value=1)
    :param rho_red : surface Reflectances Red  band [dl]
    :param rho_nir: surface Reflectances  NIR band [dl]
    :return: turbidity in FNU
    """

    limit_inf, limit_sup = 0.05, 0.07
    a_low, c_low = 228.1, 0.1641
    a_high, c_high = 3078.9, 0.2112

    t_low = nechad(Red, a_low, c_low)
    t_high = nechad(Nir2, a_high, c_high)
    w = (Red - limit_inf) / (limit_sup - limit_inf)
    t_mixing = (1 - w) * t_low + w * t_high

    t_low[Red >= limit_sup] = t_high[Red >= limit_sup]
    t_low[(Red >= limit_inf) & (Red < limit_sup)] = t_mixing[(Red >= limit_inf) & (Red < limit_sup)]
    t_low[t_low > 4000] = 0
    return t_low

# CDOM Brezonik et al. 2005
def cdom_brezonik(Blue, RedEdg2):


    cdom = np.exp(1.872 - 0.830 * np.log(Blue/RedEdg2))

    return cdom


functions = {
    'SPM_Nechad': {
        'function': nechad,
        'units': 'mg/l'
    },
    
    'TURB_Dogliotti': {
        'function': turb_dogliotti,
        'units': 'FNU'
    },

    'CHL_Gitelson': {
        'function': chl_gitelson2,
        'units': 'mg/m³'
    },

    'CDOM_Brezonik': {
        'function': cdom_brezonik,
        'units': '',
    }
}
