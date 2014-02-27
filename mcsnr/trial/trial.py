import numpy as np
import scipy.ndimage as nd
import matplotlib.pylab as plt
from astropy.table import Table

import geminiutil.gmos as gmos
import geminiutil.gmos.basic.wavecal as wavecal
import geminiutil.gmos.basic.estimate_disp as estimate_disp
import extract_psf.extract as extract


plt.ion()


do_prep = False
do_extract = False
do_calibrate = False
do_plot = False

proj = gmos.GMOSMOSProject('sqlite:///gmos.db3', work_dir='.')
sci_set = proj.science_sets[0]
# sci_set.science.prepare_to_database()
# sci_set.calculate_slice_geometries_to_database()
priority1 = [slc for slc in sci_set.slices if slc.priority == 1]
specs = []
for slc in priority1:
    specs.append(slc.extract_point_source())



#get arc data, need to calibrate
# priority1[0].get_prepared_arc_data()[0].astype('<f4')
# specpos_x = priority1[1].science_set.science.mask.table[priority1[1].list_id]['specpos_x']
if do_prep:
    do_prepare = gmos.GMOSPrepareScienceFrame(bias_subslice=[slice(None),
                                                             slice(1,11)])
    do_prepare(raw_sci)

if do_extract:
    prep_sci_fits = raw_sci.raw_fits.prepared_fits.fits.fits_data

    instrument_setup = raw_sci.raw_fits.instrument_setup
    x_binning = instrument_setup.x_binning
    y_binning = instrument_setup.y_binning
    grating_wavelength = instrument_setup.grating_central_wavelength_value
    anamorphic_factor = instrument_setup.anamorphic_factor

    # this depends on how much of data frame was cut off
    # for 6-exp, from prepare.multiext_header_value(image, 'DATASEC')
    x_offset = 0
    ndisp = prep_sci_fits[1].data.shape[-1]
    x = x_offset + x_binning*np.arange(0.5, ndisp) - 0.5

    # for 6-amp, from prepare.multiext_header_value(image, 'RDNOISE')
    read_noise = np.array([amp.header['RDNOISE']
                           for amp in prep_sci_fits[1:]]).reshape(-1,1,1)
    ff_noise = 0.03
    model_errors = 1
    # get FWHM from seeing?
    psf = extract.PSF(form=0, guesses=[[0., 0.], 8., (2.5, '@')])

    # corresponding day arc, not yet automatic
    lines = wavecal.read('as20121011s0230.hdf5', path='lines')
    dayarc = Table.read('as20121011s0230.hdf5', path='arc')
    filt = np.array([0.15,0.85,0.0,0.06,.88,0.06,0.,0.85,0.15])
    fs = nd.convolve1d(dayarc['f'], filt, axis=0, mode='nearest')
    model_wide_slit_arc = Table([dayarc['w'].flatten(), fs.flatten()],
                                names=['w','f'])
    model_wide_slit_arc.sort('w')

    # expected position where the grating wavelength will fall for a slit
    # at specpos_x=0
    xref0 = 1024.  # lines.meta['par'][0]
    refwave_arc = lines.meta['refwave']
    refwave_sci = grating_wavelength * 10.
    # later may want to switch to
    xref0 = lines.meta['par'][0] + 15. - (refwave_sci - refwave_arc) / \
        lines.meta['par'][3]

    # loop over all slits, extracting priority=1 spectra
    for slitno, sliceinfo in enumerate(raw_sci.slices):
        if sliceinfo.priority != 1:
            continue

        sci_slice = slice(*(np.ceil(edge).astype(int)
                            for edge in (sliceinfo.lower_edge,
                                         sliceinfo.upper_edge)))

        scidata = np.array([prep_sci_fits[chip].data[sci_slice, :]
                            for chip in range(1,4)])

        tracepos = np.array([scidata.shape[1]/2. -
                             raw_sci.raw_fits.mask.table[slitno]['specpos_y'] /
                             y_binning])
        e_source = scidata
        for i in range(model_errors+1):
            # if we want to use model errors, first get approximate model
            # (iteratively, if model_errors>1)
            error_estimate = extract.ccd_noise(e_source, read_noise, ff_noise)
            (out,eout,back,chi2,test,
             ntdiscard,ntbadl,ntbadh,
             nproblems) = extract.extract(scidata, tracepos, psf,
                                          e=error_estimate, skypol=0,
                                          ibadlimit=2, squeeze_dims=False,
                                          itesttype=101 if i < model_errors
                                          else 103)
            # set error source to test frame
            e_source = test

        # make table of extracted spectra
        # Note: out.shape=(nstar,norder,nwav) but need
        #       (nwav,nstar,norder) for nstar>1 -> transpose(2,0,1)
        #       (nwav,norder) for nstar=1 -> transpose(1,0)
        scitab = Table([np.array(x, dtype=np.float32)] +
                       [a.transpose(2,0,1) for a in [out, eout, back, chi2]],
                       names=['x', 'f', 'e', 's', 'chi2'])

        # fits0_header = prep_sci_fits[0].header
        # for hdr in fits0_header:
        #     if hdr not in ('SIMPLE', 'BITPIX', 'EXTEND', ''):
        #         scitab.meta[hdr.lower()] = fits0_header[hdr]
        scitab.meta['nstars'] = len(tracepos)
        scitab.meta['tracepos'] = tracepos
        scitab.meta['psfform'] = psf.form
        scitab.meta['psfnpoldisp'] = psf.npoldisp
        scitab.meta['psfpar'] = psf.par
        scitab.meta['psferr'] = psf.err
        scitab.meta['psfchi2'] = psf.chi2
        scitab.meta['ndiscard'] = ntdiscard
        scitab.meta['nbad'] = [ntbadl, ntbadh]
        scitab.meta['nproblems'] = nproblems

        plt.ion()
        plt.clf()
        plt.plot(scitab['x'], scitab['f'][:,0,:])
        plt.draw()

if __name__ != '__main__':
        # get corresponding arc
        prep_arc_fits = raw_sci.mask_arc.prepared_fits.fits.fits_data
        arcdata = np.array([prep_arc_fits[chip].data[sci_slice, :]
                            for chip in range(1,4)])

        ix = np.round(psf.offsets()+tracepos[0]).astype(np.int)
        if len(tracepos) == 1:
            ix = ix[np.newaxis,:]
        (_, arcchi2, arctest, ntbadl, ntbadh) = extract.fitsky(arcdata)
        scitab['a'] = arctest[(0,1,2),ix[:,(0,1,2)],:].transpose(2,0,1)

        # plt.plot(dayarc['w'], dayarc['f']-5000.)
        # plt.plot(model_wide_slit_arc['w'], model_wide_slit_arc['f']-3500.)
        # scitab['w'] = lines.fit(np.arange(3)[np.newaxis,:],
        #                         scitab['x'][:,np.newaxis])[:,np.newaxis,:]
        # plt.plot(scitab['w'][:,0,:], scitab['a'][:,0,:])

        # **** AM HERE ***
        specpos_x = raw_sci.raw_fits.mask.table[slitno]['specpos_x']
        # print slitno, specpos_x
        #   specpos_x  xrefbest
        # 18:  -56. -> 624.
        # 19:  124. -> 466.
        # 20: -119  -> 679.
        # 21: -135  -> 693.
        # 22:  101  -> 487.
        # xrefguess = lines.meta['par'][0] + 236./lines.meta['par'][3]
        xrefguess = xref0 - specpos_x / anamorphic_factor

        xrefbest, dwbest, shift = estimate_disp.estimate_disp(
            Table([scitab['a'][:,0,1], scitab['x']], names=('f','x')),
            refwave_arc,
            xrefguess+np.arange(-30,30),
            lines.meta['par'][3]*np.linspace(0.97,1.03,31),
            comparison=model_wide_slit_arc, full=True)

        lines.meta['guess'] = np.hstack([xrefbest, lines.meta['par'][1:]])
        scitab['w'] = lines.estimate(np.arange(3)[np.newaxis,:],
                                     scitab['x'][:,np.newaxis])[:,np.newaxis,:]

        # plt.plot(scitab['w'][:,0,:], scitab['a'][:,0,:] + 5000.)
        # scitab.write('scitab.hdf5', path='spectrum', overwrite=True)
        plt.ion()
        plt.plot(scitab['w'][:,0,:], scitab['f'][:,0,:])
        plt.draw()
