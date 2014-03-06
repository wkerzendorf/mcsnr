import numpy as np
import scipy.ndimage as nd
import astropy.constants.si as si
import astropy.units as u
from astropy.table import Table, Column
import warnings

# unit stuff insufficiently well developped to use beyond this
def vbyc(v):
    if not hasattr(v,'unit'): v *= u.km/u.s
    return (v.si/si.c).value

def doppler(v):
    return 1.+vbyc(v)

def doppler_undo(v):
    return 1./doppler(v)

doppler_logw = vbyc

def fittable(obs, model, *args, **kwargs):
    return fit(obs['w'], obs['flux'], obs['err'], model['w'], model['flux'],
               *args, **kwargs)

def fit(w, f, e, mw, mf, vgrid, npol,
        sigrange=None, vrange=None, doppler=doppler, plot=False):

    chi2 = Table([vgrid, np.zeros_like(vgrid)], names=['v','chi2'])
    chi2['v'].unit = u.km/u.s

    Vp = np.polynomial.legendre.legvander(w/w.mean()-1., npol)

    rcond = len(f)*np.finfo(f.dtype).eps
    e2 = e**2
    fbye2 = f/e2
    for i,v in enumerate(chi2['v']):
        mfi = np.interp(w, mw*doppler(v), mf)
        V = Vp*(mfi/e2)[:,np.newaxis]
        scl = np.sqrt((V*V).sum(0))
        sol, resids, rank, s = np.linalg.lstsq(V/scl, fbye2, rcond)
        sol = (sol.T/scl).T
        if rank != npol:
            msg = "The fit may be poorly conditioned"
            warnings.warn(msg)
        fit = np.dot(V, sol)*e2
        chi2['chi2'][i] = np.sum(((f-fit)/e)**2)

    chi2.meta['ndata'] = len(f)
    chi2.meta['npar'] = npol+1+1
    chi2.meta['ndof'] = chi2.meta['ndata']-chi2.meta['npar']

    if plot:
        import matplotlib.pylab as plt
        plt.scatter(chi2['v'], chi2['chi2'])

    if vrange is None and sigrange is None:
        return chi2

    vbest, verr = minchi2(chi2, vrange, sigrange, plot=plot)

    mfi = np.interp(w, mw*doppler(vbest), mf)
    V = Vp*(mfi/e2)[:,np.newaxis]
    scl = np.sqrt((V*V).sum(0))
    sol, resids, rank, s = np.linalg.lstsq(V/scl, fbye2, rcond)
    sol = (sol.T/scl).T
    if rank != npol:
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg)
    chi2.meta['wmean'] = w.mean()
    chi2.meta['continuum'] = sol
    fit = np.dot(V, sol)*e2
    return chi2, fit, mfi

def minchi2(chi2, vrange=None, sigrange=None,
            fitcol='chi2fit', plot=False):
    assert vrange is not None or sigrange is not None
    if sigrange is None: sigrange = 1e10
    if vrange is None: vrange = 1e10

    iminchi2 = chi2['chi2'].argmin()
    ndof = float(chi2.meta['ndof'])
    iok = np.where(np.logical_and(chi2['chi2'] <
                                  chi2['chi2'][iminchi2]*(1.+sigrange**2/ndof),
                                  abs(chi2['v']-chi2['v'][iminchi2]) <=
                                  vrange))

    p = np.polynomial.Polynomial.fit(chi2['v'][iok], chi2['chi2'][iok],
                                     2, domain=[])

    if plot:
        import matplotlib.pylab as plt
        plt.scatter(chi2['v'][iok], chi2['chi2'][iok], c='g')
        plt.plot(chi2['v'], p(chi2['v']))

    vbest = -p.coef[1]/2./p.coef[2]
    # normally, get sigma from where chi2 = chi2min+1, but best to scale
    # errors, so look for chi2 = chi2min*(1+1/ndof) ->
    # a verr**2 = chi2min/ndof -> verr = sqrt(chi2min/ndof/a)
    verr = np.sqrt(p(vbest)/p.coef[2]/ndof)
    chi2.meta['vbest'] = vbest
    chi2.meta['verr'] = verr
    chi2.meta['bestchi2'] = p(vbest)
    if fitcol is not None:
        if fitcol in chi2.colnames:
            chi2[fitcol] = p(chi2['v'])
        else:
            chi2.add_column(Column(data=p(chi2['v']), name=fitcol))

    return vbest, verr

def observe(model, wgrid, slit, seeing, overresolve, offset=0.):
    """Convolve a model with a seeing profile, truncated by a slit, and pixelate

    Parameters
    ----------
    model: Table (or dict-like)
       Holding wavelengths and fluxes in columns 'w', 'flux'
    wgrid: array
       Wavelength grid to interpolate model on
    slit: float
       Size of the slit in wavelength units
    seeing: float
       FWHM of the seeing disk in wavelength units
    overresolve: int
       Factor by which detector pixels are overresolved in the wavelength grid
    offset: float, optional
       Offset of the star in the slit in wavelength units (default 0.)

    Returns
    -------
    Convolved model: Table
       Holding wavelength grid and interpolated, convolved fluxes
       in columns 'w', 'flux'
    """
    # make filter
    wgridres = np.min(np.abs(np.diff(wgrid)))
    filthalfsize = np.round(slit/2./wgridres)
    filtgrid = np.arange(-filthalfsize,filthalfsize+1)*wgridres
    # sigma ~ seeing-fwhm/sqrt(8*ln(2.))
    filtsig = seeing/np.sqrt(8.*np.log(2.))
    filt = np.exp(-0.5*((filtgrid-offset)/filtsig)**2)
    filt /= filt.sum()
    # convolve with pixel width
    filtextra = int((overresolve-1)/2+0.5)
    filt = np.hstack((np.zeros(filtextra), filt, np.zeros(filtextra)))
    filt = nd.convolve1d(filt, np.ones(overresolve)/overresolve)
    mint = np.interp(wgrid, model['w'], model['flux'])
    mconv = nd.convolve1d(mint, filt)
    return Table([wgrid, mconv], names=('w','flux'), meta={'filt': filt})
