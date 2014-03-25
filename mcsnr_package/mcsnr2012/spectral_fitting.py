import warnings
import numpy as np

from astropy import units as u, constants as const

from geminiutil.gmos.alchemy.mos import MOSSpectrum



def get_model_spectrum(self, teff, logg, feh):
   if self.spectral_grid is None:
       raise AttributeError('No spectral grid associated')

   self.spectral_grid.teff = teff
   self.spectral_grid.logg = logg
   self.spectral_grid.feh = feh

   return self.spectral_grid.wave, self.spectral_grid(self.spectral_grid.wave)


def get_spectral_fit(self, teff, logg, feh, velocity=0, npol=3):

   model_wavelength, model_flux = self.get_model_spectrum(teff, logg, feh)

   # pre-calculate observed flux/error and polynomial bases
   if not hasattr(self, 'signal_to_noise'):
       self.signal_to_noise = (self.flux / self.uncertainty)
   # Vp[:,0]=1, Vp[:,1]=w, .., Vp[:,npol]=w**npol
   if not hasattr(self, '_Vp') or self._Vp.shape[-1] != npol + 1:
       self._Vp = np.polynomial.polynomial.polyvander(
           self.wavelength/self.wavelength.mean() - 1., npol)

   return self._spectral_fit(model_wavelength, model_flux, velocity)


def _spectral_fit(self, model_wavelength, model_flux, velocity):
   # Doppler shift the model grid
   velocity = u.Quantity(velocity, u.Unit('km/s'))
   shifted_wavelength = (model_wavelength.to(self.wavelength.unit) *
                         (1 + velocity/const.c))

   # interpolate Doppler-shifted model on the observed wavelengths
   interpolated_model = np.interp(self.wavelength.value,
                                  shifted_wavelength.value,
                                  model_flux)

   rcond = len(self.flux)*np.finfo(self.flux.dtype).eps
   # V[:,0]=mfi/e, Vp[:,1]=mfi/e*w, .., Vp[:,npol]=mfi/e*w**npol
   V = self._Vp * (interpolated_model/self.uncertainty)[:,np.newaxis]
   # normalizes different powers
   scl = np.sqrt((V*V).sum(0))
   sol, resids, rank, s = np.linalg.lstsq(V/scl, self.signal_to_noise,
                                          rcond)
   sol = (sol.T/scl).T
   if rank != self._Vp.shape[-1] - 1:
       msg = "The fit may be poorly conditioned"
       warnings.warn(msg)

   fit = np.dot(V, sol) * self.uncertainty
   chi2 = np.sum(((self.flux-fit)/self.uncertainty)**2)
   return chi2, fit, interpolated_model


MOSSpectrum.get_spectral_fit = get_spectral_fit
MOSSpectrum._spectral_fit = _spectral_fit
MOSSpectrum.get_model_spectrum= get_model_spectrum