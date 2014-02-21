from __future__ import print_function, division

import numpy as np

from astropy.table import Table
import astropy.units as u

import matplotlib.pylab as plt

from geminiutil.gmos import GMOSMOSProject, GMOSPrepareFrame
from geminiutil.gmos.util.longslit_arc import GMOSLongslitArc

#import extract_psf.extract as extract


def line_catalog(lincat='/home/mhvk/standard/lines/CuAr.dat'):
    return Table(
        np.genfromtxt(lincat,
                      dtype=[('w','f8'), ('ion','a7'), ('strength','i4')],
                      delimiter=[9, 7, 8]))

if __name__ == '__main__':
    plt.ion()

    # os.system('rm mcsnr.db3')
    dbname = 'sqlite:///mcsnr.db3'
    proj = GMOSMOSProject(dbname)
    # proj.initialize_database()
    # proj.add_directory('/raw/mhvk/gemini/mcsnr', file_filter='S*S*.fits')
    # proj.add_directory('/raw/mhvk/gemini/mcsnr/mdf_dir')
    # proj.link_masks()
    # check proj.science_instrument_setups and
    # proj.longslit_arcs_instrument_setups and
    # find the best matches (science:longslit dict)
    # proj.link_science_sets({10:16, 11:16, 12:15, 13:15})

    # also have 1.5" arcs, which are easier to use as wide-slit
    # ignore for now
    # dayarc_longslit = [fil for fil in proj.daycal
    #                    if(fil.object.name == 'cuar' and
    #                       'arcsec' in fil.fits.header['maskname'])]
    # daycalfits = [fil for fil in proj.daycal
    #               if(fil.object.name == 'cuar' and
    #                  fil.fits.header['maskname'] == '0.5arcsec')]
    daycalfits = set([s.longslit_arc for s in proj.science_sets])
    bluedaycal = [fil for fil in daycalfits
                  if fil.instrument_setup.grating.name[:1] == 'B']
    reddaycal = [fil for fil in daycalfits
                 if fil.instrument_setup.grating.name[:1] == 'R'][0]
    # DARN: first & last of reddaycal taken at same wavelength

    lincat = line_catalog()

    doblue = True
    dored = False
    doarcplot = True

    prepdayarc = GMOSPrepareFrame(bias_subslice=[slice(None), slice(1,11)],
                                  data_subslice=[slice(1150,1250),slice(-1)])
    if doblue:
        calbluedayarc = GMOSLongslitArc(
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
        calbluedayarc = GMOSLongslitArc(
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
