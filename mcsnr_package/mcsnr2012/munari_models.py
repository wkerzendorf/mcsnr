import numpy as np
from scipy import optimize
from collections import OrderedDict
import warnings
from astropy.constants import si
from astropy import units as u
import pylab
from pyspec import oned

def vbyc(v):
    if not hasattr(v,'unit'): v *= u.km/u.s
    return (v.si/si.c).value

def doppler(v):
    return 1.+vbyc(v)

def doppler_undo(v):
    return 1./doppler(v)

doppler_logw = vbyc


class SpectrumModel(object):
    def __init__(self, munari_grid, spectrum):
        self.grid = munari_grid
        self.spectrum = spectrum

        self.params = OrderedDict([('p0', 1e-10), ('p1',0), ('p2',0), ('p3', 0), ('p4', 0), ('offset', 0), ('teff',5780),
                                   ('logg',4.4), ('feh',0.0), ('v_rad',0.0), ('v_rot', 2.0)])
        self.bounds = ((None, None), (None, None), (None, None), (None, None), (None, None), (None, None), (4000, 40000),
                       (1, 4.5), (-3, 0.), (-500, 500), (2, 200))
        self.fit_params = np.array(self.params.values())
        self.fixed = OrderedDict([(key, False) for key in self.params.keys()])

        self.param_no = len(self.params)
        self.flux_mask = np.ones_like(self.spectrum.wave).astype(bool)
        self.flux_mask = self.spectrum.mask
        self.current_polyval = np.polynomial.legendre.legval

    def __call__(self, p0, p1, p2, p3, p4, offset, teff, logg, feh, v_rad, v_rot):
        new_model = self.grid.interpolate_spectrum(teff=teff, logg=logg, feh=feh)
        new_model = new_model.shift_velocity(v_rad)
        new_model = new_model.convolve_rotation(abs(v_rot))

        poly_param = (p0, p1, p2, p3, p4)
        continuum_poly = self.current_polyval(new_model.wave, poly_param)
        new_model.flux *= continuum_poly
        new_model.flux += offset
        new_model = new_model.interpolate(self.spectrum.wave)
        return new_model


    def eval_current(self):
        for i, value in enumerate(self.params.values()):
            self.fit_params[i] = value

        return self(*self.fit_params)

    def chi2(self, param):
        self.fit_params[~self.mask] = np.array(param)
        #print param
        #print self.fit_params
        return np.sum(((self(*self.fit_params) - self.spectrum).flux**2)[self.flux_mask])

    def fit(self, method='Nelder-Mead'):
        self.mask = np.array(self.fixed.values())
        x0 = np.array(self.params.values())[~self.mask]
        self.fit_params = np.array(self.params.values())
        print "beginning with %s" % x0
        print "my current _params", self.params
        fit_result = optimize.minimize(self.chi2, x0, method=method, bounds=self.bounds)
        for fit_param, param_name, cur_mask in zip(self.fit_params, self.params.keys(), self.mask):
            if not cur_mask:
                self.params[param_name] = fit_param
        self.fit_result = fit_result
        return fit_result

    def do_continuum_fit(self):
        self.fixed['teff']=True
        self.fixed['logg']=True
        self.fixed['feh']=True
        self.fixed['v_rad']=True
        self.fixed['v_rot']=True
        self.fixed['offset']=True

        self.fit(method='Nelder-Mead')
        self.fit(method='Powell')



    def do_complete_fit(self, fix_vrad=False):
        self.fixed['teff']=True
        self.fixed['logg']=True
        self.fixed['feh']=True
        self.fixed['v_rad']=True
        self.fixed['v_rot']=True

        self.fit(method='Powell')
        if not fix_vrad:
            self.fixed['v_rad']=False

            self.fit(method='Powell')

        self.fixed['p0']=True
        self.fixed['p1']=True
        self.fixed['p2']=True
        self.fixed['teff']=False
        self.fixed['logg']=False
        #extr_model.fixed['feh']=False
        #extr_model.fixed['v_rot']=False
        self.fit(method='Powell')

class SpectrumModel2(object):
    def __init__(self, munari_grid, spectrum, error=None, npol=5, warn_conditioned=False):
        self.grid = munari_grid
        self.spectrum = spectrum
        self.warn_conditioned = warn_conditioned

        self.mask_regions = self.spectrum.wave[np.diff(self.spectrum.mask)]
        if error is None:
            self.error2 = np.ones_like(self.spectrum.flux)
            #self.error2 = np.sqrt(np.abs(self.spectrum.flux))**2
        else:
            self.error2 = error**2

        self.fluxbye2 = self.spectrum.flux / self.error2
        self.params = OrderedDict([('offset', 0), ('teff',5780),
                                   ('logg',4.4), ('feh',0.0), ('v_rad',0.0)])
        self.bounds = ((None, None), (4000, 40000),
                       (1, 4.7), (-3, 0.), (-500, 500))

        self.fit_params = np.array(self.params.values())
        self.fixed = OrderedDict([(key, False) for key in self.params.keys()])
        self.mask = np.array(self.fixed.values())

        self.npol = npol
        self.param_no = len(self.params)
        self.flux_mask = np.ones_like(self.spectrum.wave).astype(bool)
        self.flux_mask = self.spectrum.mask
        self.Vandermonde = np.polynomial.legendre.legvander(self.spectrum.wave/self.spectrum.wave.mean() - 1.0, npol)
        self.rcond = len(self.spectrum.flux)*np.finfo(self.spectrum.flux.dtype).eps
        self.fitspec = None



    def __call__(self, offset, teff, logg, feh, v_rad):
        new_model = self.grid.interpolate_spectrum(teff=teff, logg=logg, feh=feh)
        new_model = new_model.shift_velocity(v_rad, interp=False)
        if abs(new_model.flux[0] + 1) < 1e-10 :
            print 'outside of grid'
            new_model.flux *= 1e10
            return new_model.interpolate(self.spectrum.wave)

        new_model_flux = np.interp(self.spectrum.wave, new_model.wave, new_model.flux) + offset

        V = self.Vandermonde * (new_model_flux/ self.error2)[:,np.newaxis]

        scl = np.sqrt((V * V).sum(0))

        sol, resids, rank, s = np.linalg.lstsq(V/scl, self.fluxbye2, self.rcond)
        sol = (sol.T/scl).T
        if (rank != self.npol) and self.warn_conditioned:
            msg = "The fit may be poorly conditioned"
            warnings.warn(msg)
        fit = np.dot(V, sol) * self.error2
        fitspec = oned.onedspec(self.spectrum.wave, fit, mode='waveflux')
        self.fitspec = fitspec
        return fitspec




    def eval_current(self):
        for i, value in enumerate(self.params.values()):
            self.fit_params[i] = value
        self.fitspec = self(*self.fit_params)
        return self.fitspec

    def chi2(self, param):
        self.fit_params[~self.mask] = np.array(param)
        #print param
        #print self.fit_params

        return np.sum(((self(*self.fit_params).flux - self.spectrum.flux)**2)[self.flux_mask])

    def plot(self, fig):
        if self.fitspec is None:
            raise ValueError('fitspec None - run fit first')
        fig.clf()

        ax = fig.add_subplot(111)
        ax.plot(self.spectrum.wave, self.spectrum.flux)
        for mask_start, mask_end in zip(self.mask_regions[:-1][::2], self.mask_regions[1:][::2]):
            ax.axvspan(mask_start, mask_end, color='red', alpha=0.3)
        ax.plot(self.fitspec.wave, self.fitspec.flux)
        ax.set_title('teff=%.2f logg=%.2f feh=%.2f vrad=%g' % \
                     (self.params['teff'], self.params['logg'], self.params['feh'], self.params['v_rad']))


    def grid_fit(self):
        pass

    def fit(self, method='Nelder-Mead'):
        self.mask = np.array(self.fixed.values())
        x0 = np.array(self.params.values())[~self.mask]
        self.fit_params = np.array(self.params.values())
        print "beginning with %s" % x0
        print "my current _params", self.params
        current_bounds = [bound for (bound, fixed) in zip(self.bounds, self.fixed.values()) if not fixed]
        fit_result = optimize.minimize(self.chi2, x0, method=method, bounds=current_bounds)
        for fit_param, param_name, cur_mask in zip(self.fit_params, self.params.keys(), self.mask):
            if not cur_mask:
                self.params[param_name] = fit_param
        self.fit_result = fit_result
        return fit_result


    def grid_fit(self, v_rads, teffs, loggs, v_rad_sample=150., teff_sample=500, logg_sample=1., plot_fit=False):
        print "First doing the vrads"
        v_rads_chi2 = np.zeros_like(v_rads)
        for i, v_rad in enumerate(v_rads):
            self.params['v_rad'] = v_rad
            v_rads_chi2[i] = self.chi2(self.params.values())

        v_chi2_min = v_rads[np.argmin(v_rads_chi2)]
        v_mask = (v_rads > (v_chi2_min - v_rad_sample)) & (v_rads < (v_chi2_min + v_rad_sample))
        p = np.polynomial.Polynomial.fit(v_rads[v_mask], v_rads_chi2[v_mask], 2, domain=[])
        v_best = -p.coef[1]/2./p.coef[2]

        self.params['v_rad'] = v_best

        if plot_fit:
            pylab.figure(1)
            pylab.clf()
            pylab.plot(v_rads, v_rads_chi2, label='chi2')
            pylab.plot(v_rads, p(v_rads), label='fit')
            pylab.axvline(v_best)
            pylab.title('v_best=%g' % v_best)
            pylab.legend()

        print "Fitting teffs, logg grid"

        teff_logg_grid = np.zeros((len(loggs), len(teffs)))
        for i, teff, in enumerate(teffs):
            for j, logg in enumerate(loggs):
                self.params['teff'] = teff
                self.params['logg'] = logg
                teff_logg_grid[j, i] = self.chi2(self.params.values())
        teff_logg_grid = np.ma.MaskedArray(teff_logg_grid, mask=(teff_logg_grid>1e13) )


        teff_chi2 = teff_logg_grid.min(axis=0)
        logg_chi2 = teff_logg_grid.min(axis=1)

        teff_chi2_min = teffs[np.argmin(teff_chi2)]
        logg_chi2_min = loggs[np.argmin(logg_chi2)]

        teff_mask = (teffs > (teff_chi2_min - teff_sample)) & (teffs < (teff_chi2_min + teff_sample))
        logg_mask = (loggs > (logg_chi2_min - logg_sample)) & (loggs < (logg_chi2_min + logg_sample))

        teff_p = np.polynomial.Polynomial.fit(teffs[teff_mask], teff_chi2[teff_mask], 2, domain=[])
        logg_p = np.polynomial.Polynomial.fit(loggs[logg_mask], logg_chi2[logg_mask], 2, domain=[])

        teff_best = -teff_p.coef[1]/2./teff_p.coef[2]
        logg_best = -logg_p.coef[1]/2./logg_p.coef[2]

        if plot_fit:
            pylab.figure(2)
            pylab.clf()
            pylab.imshow(teff_logg_grid, extent=(teffs.min(), teffs.max(), loggs.min(), loggs.max()), aspect='auto', interpolation='nearest')
            pylab.plot([teff_best], [logg_best], 'rx', mew=3, ms=25)
            pylab.xlabel('Teff')
            pylab.ylabel('log(g)')
            pylab.colorbar()



        return v_rads_chi2, v_best, teff_best, logg_best, teff_logg_grid


    def do_continuum_fit(self):
        self.fixed['teff']=True
        self.fixed['logg']=True
        self.fixed['feh']=True
        self.fixed['v_rad']=True
        self.fixed['v_rot']=True
        self.fixed['offset']=True

        self.fit(method='Nelder-Mead')
        self.fit(method='Powell')



    def do_complete_fit(self, fix_vrad=False):
        self.fixed['teff']=True
        self.fixed['logg']=True
        self.fixed['feh']=True
        self.fixed['v_rad']=True
        self.fixed['v_rot']=True

        self.fit(method='Powell')
        if not fix_vrad:
            self.fixed['v_rad']=False

            self.fit(method='Powell')

        self.fixed['p0']=True
        self.fixed['p1']=True
        self.fixed['p2']=True
        self.fixed['teff']=False
        self.fixed['logg']=False
        #extr_model.fixed['feh']=False
        #extr_model.fixed['v_rot']=False
        self.fit(method='Powell')



def gaussian1d(x, peak, x0, sigma_x, offset):
    y = abs(offset) + abs(peak) * np.exp(-(x - x0)**2 / (2**.5 * sigma_x))
    return y