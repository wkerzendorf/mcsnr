import numpy as np
import scipy.ndimage as nd
import matplotlib.pylab as plt
import astropy.units as u
from astropy.table import Table
from scipy.optimize import leastsq
import numpy.polynomial.polynomial as poly

from geminiutil.gmos import gmos_alchemy as ga
import geminiutil.gmos as gmos
import geminiutil.gmos.util.wavecal as wavecal
import geminiutil.gmos.util.estimate_disp as estimate_disp
import extract_psf.extract as extract


plt.ion()


do_prep = False
do_extract = False
do_calibrate = False
do_plot = False


def selpar(guess):
    # return np.hstack([par[0], par[3:]])
    return np.hstack([guess[0], guess[1][3:]])


def allpar(par):
    return par[0], np.hstack([guess[1][:3], par[1:]])


def residual(par, model):
    newlines.meta['refwave'], newlines.meta['par'] = allpar(par)
    wguess = newlines.fit(np.arange(spec['a'].shape[-1]).reshape(1,1,-1),
                          spec['x'].reshape(-1,1,1))
    mguess = np.interp(wguess, model['w'], model['f'])
    pfit, extra = poly.polyfit(mguess.flatten(), spec['a'].flatten(),
                               1, full=True)
    print('par={}, chi2={}'.format(par, extra[0]))
    aguess = pfit[0] + pfit[1]*mguess
    return (spec['a'] - aguess).flatten()


proj = gmos.GMOSMOSProject('sqlite:///gmos.db3', work_dir='.')
# proj.metadata.create_all()

# daycals = {'B': {}, 'R': {}}
# daycalfits = proj.session.query(ga.GMOSLongSlitArc).all()
# daycalfits = [fil for fil in daycalfits
#               if fil.raw.mask.name.startswith('0.5')]
# for daycal in daycalfits:
#     daycals[daycal.raw.instrument_setup.grating.name[:1]][
#         daycal.raw.instrument_setup.grating_slit_wavelength_value
#     ] = daycal.wave_cal.fname

for sci_set in proj.science_sets:
    try:
        sci_set.science.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    try:
        sci_set.calculate_slice_geometries_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    try:
        sci_set.mask_arc.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    sci_set.extract_spectra_and_calibrate()


# priority1 = [slc for slc in sci_set.slices
#              if any(src.priority == 1 for src in slc.mos_point_sources)]
# lower_edge: 1128.50433357,1094.51992124,1244.68764759,1213.14774821,1171.09892997
# specpos_x:  -56.38671875,124.385742188,-119.149658203,-134.937255859,100.985839844
# 0,1: sim Y, diff X: 0 like 2,3; 1 OKish
# 1,2 diff Y,X; 1 OKish, 2 bad
# 2,3: sim X,Y -> similarly off
# 3,4: sim Y, diff X -> 3 off, 4 OKish
# adding dw -> everything worse
# using refwave instead of xref -> much more consistent.  With dw, overshoots
# priority1 = priority1[3:5]
# spectra = []

if do_calibrate:
    instrument_setup = sci_set.science.instrument_setup
    grating_wavelength = instrument_setup.grating_central_wavelength_value
    # anamorphic_factor = instrument_setup.anamorphic_factor

    # corresponding day arc, not yet automatic
    closest = min(daycals[instrument_setup.grating.name[:1]].keys(),
                  key=lambda x: abs(grating_wavelength-x))
    dayhdf5 = daycals[instrument_setup.grating.name[:1]][closest]
    lines = wavecal.read(dayhdf5, path='lines')
    dayarc = Table.read(dayhdf5, path='arc')
    filt = np.array([0.15,0.85,0.0,0.06,.88,0.06,0.,0.85,0.15])
    fs = nd.convolve1d(dayarc['f'], filt, axis=0, mode='nearest')
    model_wide_slit_arc = Table([dayarc['w'].flatten(), fs.flatten()],
                                names=['w','f'])
    model_wide_slit_arc.sort('w')
    model_wide_slit_arc.meta.update(lines.meta)

    # # expected position where the grating wavelength will fall for a slit
    # # at specpos_x=0
    # refwave_arc = model_wide_slit_arc.meta['refwave'] * u.AA
    # # later may want to switch to
    # par0 = model_wide_slit_arc.meta['par'].copy()
    # xref0 = ((par0[0] * u.pix -
    #           (grating_wavelength - refwave_arc) /
    #           (par0[3]) * u.AA / u.pix)).to(u.pix).value

    for slc in priority1:
        spectrum = slc.extract_and_calibrate_spectrum(model_wide_slit_arc)
        plt.plot(spectrum['wave'][:,0,:], spectrum['source'][:,0,:])
        spectra.append(spectrum)
        plt.draw()
