from __future__ import print_function, division

import numpy as np
import astropy.io.fits as fits
from astropy.table import Table
import astropy.units as u

import matplotlib.pylab as plt

from geminiutil.gmos import GMOSMOSProject, GMOSPrepare, GMOSDayArc

#import extract_psf.extract as extract

import helpers

if __name__ == '__main__':
    plt.ion()

    # os.system('rm mcsnr.db3')
    dbname = 'sqlite:///mcsnr.db3'
    proj = GMOSMOSProject(dbname)
    # proj.initialize_database()
    # proj.add_directory('/raw/mhvk/gemini/mcsnr', file_filter='S*S*.fits')
    # proj.add_directory('/raw/mhvk/gemini/mcsnr/mdf_dir')
    # proj.link_masks()
    # proj.link_science_frames()

    # also have 1.5" arcs, which are easier to use as wide-slit
    # ignore for now
    # dayarc_longslit = [fil for fil in proj.daycal
    #                    if(fil.object.name == 'cuar' and
    #                       'arcsec' in fil.fits.header['maskname'])]
    daycalfits = [fil for fil in proj.daycal
                  if(fil.object.name == 'cuar' and
                     fil.fits.header['maskname'] == '0.5arcsec')]
    bluedaycal = [fil for fil in daycalfits
                  if fil.instrument_setup.grating.name[:1] == 'B']
    reddaycal = [fil for fil in daycalfits
                 if fil.instrument_setup.grating.name[:1] == 'R']
    # DARN: first & last of reddaycal taken at same wavelength

    lincat = helpers.line_catalog()

    doblue = False
    dored = False
    doarcplot = True

    prepdayarc = GMOSPrepare(bias_subslice=[slice(None), slice(1,11)],
                             data_subslice=[slice(1150,1250),slice(-1)])
    if doblue:
        calbluedayarc = GMOSDayArc(
            line_catalog=lincat['w'] * u.Angstrom,
            min_curvature=[5.,3.,2.],
            minlist1=[(3.,1e3), (1.,3e2), (0.26,0.), (0.1,0.)],
            minlist2=[(2.,0.), (0.26,0.), (0.1,0.)])

        # first one works, second one not, third one does
        for fil in bluedaycal:
            prepdayarc(fil)
            arctab, linesall, shift, fake = \
                calbluedayarc(fil.prepared_fits, doplot=True)

            arcname = 'arc-{}.hdf5'.format(fil.fits.fname.split('.')[0])
            arctab.write(arcname, path='arc', overwrite=True)
            linesall.write(arcname, path='lines', append=True)

    if dored:
        calbluedayarc = GMOSDayArc(
            line_catalog=lincat['w'] * u.Angstrom,
            min_curvature=[1.,1.,0.5],
            minlist1=[(3.,0.), (1.,0.), (0.26,0.)],
            minlist2=[(2.,0.), (1.,0.), (0.26,0.)])

        for fil in reddaycal:
            prepdayarc(fil)
            arctab, linesall, shift, fake = \
                calbluedayarc(fil.prepared_fits, doplot=True)

            arcname = 'arc-{}.hdf5'.format(fil.fits.fname.split('.')[0])
            arctab.write(arcname, path='arc', overwrite=True)
            linesall.write(arcname, path='lines', append=True)

    if doarcplot:
        for fil in daycalfits:
            arcname = 'arc-{}.hdf5'.format(fil.fits.fname.split('.')[0])
            arctab = Table.read(arcname, path='arc')
            plt.plot(arctab['w'], arctab['f'])
