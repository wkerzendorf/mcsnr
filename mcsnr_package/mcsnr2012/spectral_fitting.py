import numpy as np

from astropy import units as u, constants as const

def get_spectral_fit(self, teff, logg, feh, velocity=0):

    velocity = u.Quantity(velocity, u.Unit('km/s'))
    if self.spectral_grid is None:
        raise AttributeError('No spectral grid associated')

    self.spectral_grid.teff = teff
    self.spectral_grid.logg = logg
    self.spectral_grid.feh = feh

    grid_flux = self.spectral_grid(self.wavelength.value)
    doppler_shift = 1 + (velocity/const.c)

    new_wavelength = self.spectral_grid.wave.value * doppler_shift.value


    return np.interp(self.wavelength.value, new_wavelength, grid_flux)
