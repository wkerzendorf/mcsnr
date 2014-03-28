import warnings
import numpy as np

from astropy import units as u, constants as const
from astropy.table import Table

from geminiutil.gmos.alchemy.mos import MOSSpectrum


def get_model_spectrum(self, teff, logg, feh):
    if self.spectral_grid is None:
        raise AttributeError('No spectral grid associated')

    self.spectral_grid.teff = teff
    self.spectral_grid.logg = logg
    self.spectral_grid.feh = feh

    return self.spectral_grid.wave, self.spectral_grid(self.spectral_grid.wave)

def get_spectral_fit(self, teff0, logg0, feh0, vrad0, vrot0, npol=5):
    spectral_model_fit = SimpleStellarParametersFit(self.spectral_grid, self)

    #optimize.minimize(spectral_model_fit, (teff0, logg0, feh0, vrad0, vrot0))




class SimpleStellarParametersFit(object):

    def __init__(self, model, spectrum, npol=5):
        self.model = model
        self.spectrum = spectrum
        self.spectrum.signal_to_noise = (self.spectrum.flux /
                                         self.spectrum.uncertainty)
        self.spectrum._Vp = np.polynomial.polynomial.polyvander(
            self.spectrum.wavelength/self.spectrum.wavelength.mean() - 1.,
            npol)
        self.spectrum._normalize = _normalize

    def __call__(self, teff, logg, feh, vrad, vrot):
        model_spec = self.model.eval(self.spectrum, teff=teff, logg=logg,
                                     feh=feh, vrad=vrad, vrot=vrot)

        normalized_model = self.spectrum._normalize(model_spec.wavelength,
                                                    model_spec.flux)
        return ((self.spectrum.flux - normalized_model) /
                self.spectrum.uncertainty)**2


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
    return fit, chi2, interpolated_model


def _normalize(self, model_wavelength, model_flux):
    # interpolate Doppler-shifted model on the observed wavelengths
    interpolated_model = np.interp(self.wavelength.value,
                                   model_wavelength.value,
                                   model_flux.value)

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
    # chi2 = np.sum(((self.flux-fit)/self.uncertainty)**2)
    return fit


def find_velocity(self, teff, logg, feh,
                  velocities=np.linspace(-400, 400, 51) * u.km / u.s,
                  npol=3, sigrange=None, vrange=70):

    teff = min(max(teff, 4000), 35000)

    model_wavelength, model_flux = self.get_model_spectrum(teff, logg, feh)
    # pre-calculate observed flux/error and polynomial bases
    if not hasattr(self, 'signal_to_noise'):
        self.signal_to_noise = (self.flux / self.uncertainty)
    # Vp[:,0]=1, Vp[:,1]=w, .., Vp[:,npol]=w**npol
    if not hasattr(self, '_Vp') or self._Vp.shape[-1] != npol + 1:
        self._Vp = np.polynomial.polynomial.polyvander(
            self.wavelength/self.wavelength.mean() - 1., npol)

    chi2 = Table([velocities, np.zeros_like(velocities.value)],
                 names=['velocity','chi2'])

    chi2['chi2'] = np.array([
        self._spectral_fit(model_wavelength, model_flux, v)[1]
        for v in velocities])
    chi2.meta['ndata'] = len(self.flux)
    chi2.meta['npar'] = npol+1+1
    chi2.meta['ndof'] = chi2.meta['ndata']-chi2.meta['npar']

    if vrange is None and sigrange is None or len(velocities) < 3:
        ibest = chi2['chi2'].argmin()
        vbest, bestchi2 = chi2[ibest]
        chi2.meta['vbest'] = vbest
        chi2.meta['verr'] = 0.
        chi2.meta['bestchi2'] = bestchi2
    else:
        vbest, verr, bestchi2 = minchi2(chi2, vrange, sigrange)

    fit, bestchi2, interpolated_model = self._spectral_fit(model_wavelength,
                                                           model_flux,
                                                           vbest)

    return vbest, verr, fit, chi2, interpolated_model


def minchi2(chi2, vrange=None, sigrange=None, fitcol='chi2fit'):
    assert vrange is not None or sigrange is not None
    if sigrange is None:
        sigrange = 1e10
    if vrange is None:
        vrange = 1e10

    iminchi2 = chi2['chi2'].argmin()
    ndof = float(chi2.meta['ndof'])
    iok = np.where((chi2['chi2'] <
                    chi2['chi2'][iminchi2]*(1.+sigrange**2/ndof)) &
                   (abs(chi2['velocity']-chi2['velocity'][iminchi2]) <= vrange))

    p = np.polynomial.Polynomial.fit(chi2['velocity'][iok], chi2['chi2'][iok],
                                     2, domain=[])

    vbest = -p.coef[1]/2./p.coef[2]
    # normally, get sigma from where chi2 = chi2min+1, but best to scale
    # errors, so look for chi2 = chi2min*(1+1/ndof) ->
    # a verr**2 = chi2min/ndof -> verr = sqrt(chi2min/ndof/a)
    verr = np.sqrt(p(vbest)/p.coef[2]/ndof)
    chi2.meta['vbest'] = vbest
    chi2.meta['verr'] = verr
    chi2.meta['bestchi2'] = p(vbest)
    if fitcol is not None:
        chi2[fitcol] = p(chi2['velocity'])

    return chi2.meta['vbest'], chi2.meta['verr'], chi2.meta['bestchi2']


MOSSpectrum.get_spectral_fit = get_spectral_fit
MOSSpectrum._spectral_fit = _spectral_fit
MOSSpectrum.get_model_spectrum = get_model_spectrum
MOSSpectrum.find_velocity = find_velocity
