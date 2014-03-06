#casagrande et al. 2010 
#now in phot2stel_param.py for adding to the db
import os
import numpy as np

cur_path = os.path.dirname(__file__)
casagrande_file = os.path.join(cur_path, 'casagrande_coeff.dat')

rv=3.1
ebv = 0.11
extin = {'bv':1.0,
        'ub':rv * (1.56 - 1.31),
        'vi':rv * (1.0 - 0.479),
        'vj':rv * (1.0 - 0.282),
        'vk':rv * (1.0 - 0.108),
        'hk':rv * (0.176 -  0.108),
        'jk':rv * (0.282 -  0.108)}
        
casagrande_coeff = np.loadtxt(casagrande_file, usecols = (1,2,3,4,5,6), delimiter=',')
polys = {'bv': casagrande_coeff[0],
        'vi': casagrande_coeff[3],
        'vk': casagrande_coeff[6],
        'vj': casagrande_coeff[4]}
        
def col2teff(color, feh, coef):
    return 5040 / (coef[0] + coef[1]*color + coef[2]*color**2 + coef[3]*feh*color + coef[4]*feh + coef[5]*feh**2)
        
