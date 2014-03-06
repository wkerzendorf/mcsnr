from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, ForeignKey, Boolean, DateTime
from sqlalchemy import update
from sqlalchemy import sql
from sqlalchemy.orm import Session
from sqlalchemy.orm import sessionmaker, relationship, backref
import numpy as np
import ephem
import os
from dateutil import parser
from astropy.utils import misc
from pyspec import oned
import pandas as pd
from astropy.io import fits
from astropy import time
import h5py
from pywcsw import wcstools

Base = declarative_base()

cent_waves = {'blue': [430, 470], 'red':[735, 775]}
grating_names = {'blue': 'B600+_G5323', 'red':'R400+_G5325'}
reduction_step_dict = {'simple':'gs',
                       'wavecal':'tgs',
                       'skysub' : 'stgs',
                       'extracted': 'estgs'}

class MC(Base):
    __tablename__ = 'mc'

    id = Column(Integer, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    short = Column(String(100))
    name = Column(String(100))
    distance = Column(Float)
    def __init__(self, ra, dec, short,name, distance):
        self.ra = ra
        self.dec = dec
        self.short = short
        self.name = name
        self.distance = distance

#lmc = MC(80.8939, -69.7561, 'lmc', 'Large Magellanic Cloud')
#smc = MC(13.1866, -72.8286, 'smc', 'Small Magellanic Cloud')

class MCPS(Base):
    __tablename__ = 'mcps'

    id = Column(Integer, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    u = Column(Float)
    u_err = Column(Float)
    b = Column(Float)
    b_err = Column(Float)
    v = Column(Float)
    v_err = Column(Float)
    i = Column(Float)
    i_err = Column(Float)
    j = Column(Float)
    j_err = Column(Float)
    h = Column(Float)
    h_err = Column(Float)
    k = Column(Float)
    k_err = Column(Float)
    offset = Column(Float)
    flag = Column(Integer)
    mc = Column(Integer, ForeignKey('mc.id'))


    def get_ra_sex(self):
        return ephem.hours(np.deg2rad(self.ra))
    def get_dec_sex(self):
        return ephem.degrees(np.deg2rad(self.dec))

    def set_ra_sex(self, value):
        self.ra = ephem.hours(value)/ephem.degree
    def set_dec_sex(self, value):
        self.dec = ephem.degrees(value)/ephem.degree

    ra_sex = property(get_ra_sex, set_ra_sex)
    dec_sex = property(get_dec_sex, set_dec_sex)

    def __init__(self, ra, dec,
                 u, u_err,
                 b, b_err,
                 v, v_err,
                 i, i_err,
                 j, j_err,
                 h, h_err,
                 k, k_err,
                 offset,
                 flag,
                 mc
                 ):

        self.ra = ra
        self.dec = dec
        self.u, self.u_err = validate_mag(u, u_err)
        self.b, self.b_err = validate_mag(b, b_err)
        self.v, self.v_err = validate_mag(v, v_err)
        self.i, self.i_err = validate_mag(i, i_err)
        self.j, self.j_err = validate_mag(j, j_err)
        self.h, self.h_err = validate_mag(h, h_err)
        self.k, self.k_err = validate_mag(k, k_err)
        self.offset = offset
        self.flag = flag
        self.mc = mc

    def __repr__(self):
        return "U = %.2f +- %.2f B = %.2f +- %.2f V = %.2f +- %.2f I = %.2f +- %.2f J = %.2f +- %.2f H = %.2f +- %.2f K = %.2f +- %.2f" % \
            (self.u, self.u_err, self.b, self.b_err, self.v, self.v_err, self.i, self.i_err, self.j, self.j_err, self.h, self.h_err, self.k, self.k_err)


class GalexBand(Base):
    __tablename__ = 'galexband'
    id = Column(Integer, primary_key=True)
    short = Column(String(100))
    name = Column(String(100))

    def __init__(self, short, name):
        self.short = short
        self.name = name

band1 = GalexBand('nuv', 'NUV only detection')
band2 = GalexBand('fuv', 'FUV only detection')
band3 = GalexBand('both', 'detection in both NUV/FUV')

class Galex(Base):
    __tablename__ = 'galex'

    id = Column(Integer, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    nuv = Column(Float)
    nuv_err = Column(Float)
    fuv = Column(Float)
    fuv_err = Column(Float)
    fov_radius = Column(Float)
    band = Column(Integer, ForeignKey('galexband.id'))
    mc = Column(Integer, ForeignKey('mc.id'))

    def __init__(self, ra, dec,
                 nuv, nuv_err,
                 fuv, fuv_err,
                 fov_radius,
                 band, mc
                ):

        self.ra = ra
        self.dec = dec
        self.nuv, self.nuv_err = validate_mag(nuv, nuv_err)
        self.fuv, self.fuv_err = validate_mag(fuv, fuv_err)
        self.fov_radius = fov_radius
        self.band = band
        self.mc = mc

class UIT(Base):
    __tablename__ = 'uit'

    id = Column(Integer, primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    m162 = Column(Float)
    m162_err = Column(Float)
    b1  = Column(Float)
    b1_err  = Column(Float)
    b5  = Column(Float)
    b5_err  = Column(Float)
    mc = Column(Integer, ForeignKey('mc.id'))
    def __init__(self, ra, dec,
                 m162, m162_err,
                 b1, b1_err,
                 b5, b5_err, mc
    ):

        self.ra = ra
        self.dec = dec
        self.m162, self.m162_err = validate_mag(m162, m162_err)
        self.b1, self.b1_err = validate_mag(b1, b1_err)
        self.b5, self.b5_err = validate_mag(b5, b5_err)
        self.mc = mc

class SNRBand(Base):
    __tablename__ = 'snrband'
    id = Column(Integer, primary_key=True)
    short = Column(String(100))
    name = Column(String(100))

    def __init__(self, short, name):
        self.short = short
        self.name = name



radio = SNRBand('R', 'radio')
optical = SNRBand('O', 'optical')
xray = SNRBand('X', 'x-ray')

class SNRType(Base):
    __tablename__ = 'snr_types'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    description = Column(String)

    def __init__(self, name, description):
        self.name = name
        self.description = description



class SNRS(Base):
    __tablename__ = 'snrs'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    ra = Column(Float)
    dec = Column(Float)
    radius = Column(Float)
    band_id = Column(Integer, ForeignKey('snrband.id'))
    mc_id = Column(Integer, ForeignKey('mc.id'))
    master_id = Column(Integer, ForeignKey('master_snrs.id'))
    given_type_id = Column(Integer, ForeignKey('snr_types.id'))
    mc = relationship('MC')
    band = relationship('SNRBand')



    master = relationship('MasterSNRS', backref='badenes')
    given_type = relationship(SNRType)

    def __init__(self, name, ra, dec, radius, band_id, mc_id):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.band_id = band_id
        self.mc_id = mc_id

    def __repr__(self):
        if self.radius is not None:
            radius = '%.2f' % (self.radius)
        else:
            radius = 'NA'
        return "<SNR %s RA=%s Dec=%s Radius=%s\">" % \
            (self.name, ephem.hours(np.deg2rad(self.ra)), ephem.degrees(np.deg2rad(self.dec)), radius)

    def remap_fits(self, fname, out_fname=None, x_scale=1., y_scale=1.):
        if out_fname is None:
            out_fname = '%s-%s-cropped.fits' % (fname.lower().replace('.fits', ''), self.name)
        scale = wcstools.getScale(fname) #arcsec per pix
        syscmd = 'remap -j %(ra).5f %(dec).5f -y %(width).2f %(height).2f -o %(out_fname)s %(fname)s'
        syscmd = syscmd % dict(ra=self.ra, dec=self.dec, width=scale * self.radius * x_scale,
                               height=scale * self.radius * y_scale, fname=fname, out_fname=out_fname)

        print syscmd
        os.system(syscmd)


class ChandraSNRS(Base):
    __tablename__ = 'chandra_snrs'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    common_name = Column(String)
    link = Column(String)
    description = Column(String)
    ra = Column(Float) # degrees
    dec = Column(Float) # degrees
    width = Column(Float) # arcminutes
    height = Column(Float) #arcminutes
    fname = Column(String)
    path = Column(String)
    mc_id = Column(Integer, ForeignKey('mc.id'))
    master_id = Column(Integer, ForeignKey('master_snrs.id'))

    master = relationship('MasterSNRS', backref='chandra')

    def __init__(self, name, common_name, link, description, ra, dec, width, height, fname, path, mc_id, master_id=None):
        self.name = name
        self.common_name = common_name
        self.link = link
        self.description = description
        self.ra = ra
        self.dec = dec
        self.width = width
        self.heigh = height
        self.fname = fname
        self.path = path
        self.mc_id = mc_id
        self.master_id = master_id

class MCSNRorgSNRS(Base):
    __tablename__ = 'mcsnr_org_snrs'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    link = Column(String)
    ra = Column(Float) # degrees
    dec = Column(Float) # degrees
    mc_id = Column(Integer, ForeignKey('mc.id'))
    master_id = Column(Integer, ForeignKey('master_snrs.id'))

    master = relationship('MasterSNRS', backref='mcsnr_org')

    def __init__(self, name, link, ra, dec, mc_id, master_id=None):
        self.name = name
        self.link = link
        self.ra = ra
        self.dec = dec
        self.mc_id = mc_id
        self.master_id = master_id


class MasterSNRS(Base):
    __tablename__ = 'master_snrs'


    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float) # degrees
    dec = Column(Float) # degrees
    radius = Column(Float) #arcseconds
    width = Column(Float)
    height = Column(Float)
    mcels_available = Column(Boolean)
    chandra_available = Column(Boolean)
    mc_id = Column(Integer, ForeignKey('mc.id'))

    mc = relationship('MC')



    def __init__(self, name, ra, dec, radius, width, height):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.width = width
        self.height = height


class NewGeometry(Base):
    __tablename__ = 'new_geometry'
    id = Column(Integer, ForeignKey('snrs.id'), primary_key=True)
    ra = Column(Float)
    dec = Column(Float)
    radius = Column(Float)
    radius_fraction = Column(Float)
    pa = Column(Float)
    mag_thresh = Column(Float)
    mag_type = Column(String(100))

    snr = relationship('SNRS')
    def __init__(self, ra, dec, radius):
        self.ra = ra
        self.dec = dec
        self.radius = radius

class SNRGeminiTarget(Base):
    __tablename__ = 'snr_gemini_target'
    id = Column(Integer, primary_key=True)
    snr_id = Column(Integer, ForeignKey('snrs.id'))
    ra = Column(Float)
    dec = Column(Float)
    radius = Column(Float)
    radius_fraction = Column(Float)
    pa = Column(Float)
    mag_thresh = Column(Float)
    mag_type = Column(String(100))
    special_bv_center = Column(Float)
    special_bv_slope = Column(Float)
    special_bv_mag_width = Column(Float)
    special_bv_blue_thresh = Column(Float)
    special_bv_blue_depth = Column(Float)
    snr_location_xray_min = Column(Float)
    snr_location_xray_max = Column(Float)
    snr_location_xray_smooth = Column(Float)

    snr_location_radio_min = Column(Float)
    snr_location_radio_max = Column(Float)
    snr_location_radio_smooth = Column(Float)
    target_id = Column(Integer)
    target_priority = Column(Integer)
    observed = Column(Boolean)
    reduced = Column(Boolean)
    rereduce = Column(Boolean)


    snr = relationship('SNRS', uselist=False, backref=backref('gemini_target', uselist=False))


    mdf_dir = None
    reduced_dir = None

    reduction_step_dict = {'simple':'gs',
                           'wavecal':'tgs',
                           'skysub' : 'stgs',
                           'extracted': 'estgs'}

    _mdf = None


    def get_ra_sex(self):
        return ephem.hours(np.deg2rad(self.ra))
    def get_dec_sex(self):
        return ephem.degrees(np.deg2rad(self.dec))

    def set_ra_sex(self, value):
        self.ra = ephem.hours(value)/ephem.degree
    def set_dec_sex(self, value):
        self.dec = ephem.degrees(value)/ephem.degree

    ra_sex = property(get_ra_sex, set_ra_sex)
    dec_sex = property(get_dec_sex, set_dec_sex)



    def __init__(self, snr_id, ra, dec, radius, radius_fraction, pa=0, mag_thresh=18, mag_type='v',
                        special_bv_center=0.5, special_bv_slope=0.2, special_bv_mag_width=0.2,
                        special_bv_blue_thresh=-0.2, special_bv_blue_depth=19.5,
                        target_id=0, target_priority=0):
        self.snr_id = snr_id
        self.ra = ra
        self.dec = dec
        self.radius = radius
        self.radius_fraction = radius_fraction
        self.pa = pa,
        self.mag_thresh = mag_thresh
        self.mag_type = mag_type
        self.special_bv_blue_depth = special_bv_blue_depth
        self.special_bv_blue_thresh = special_bv_blue_thresh
        self.special_bv_center = special_bv_center
        self.special_bv_slope = special_bv_slope
        self.special_bv_mag_width = special_bv_mag_width
        self.target_id = target_id
        self.target_priority = target_priority

    def __repr__(self):
        return "<Gemini SNR #%d (%s) Target RA=%s Dec=%s Radius=%.2f\" Inner Radius Fraction=%.4f PA=%.2f>" % \
            (self.snr_id, self.snr.name, self.ra_sex, self.dec_sex,
                self.radius, self.radius_fraction, self.pa)


    @misc.lazyproperty
    def session(self):
        return Session.object_session(self)

    @property
    def science(self):
        return self.session.query(RawFITSFiles).filter_by(obs_class='science', snr_id=self.snr_id).all()

    @property
    def candidates(self):
        return self.session.query(Candidate).filter_by(snr_id=self.snr_id, priority=1).all()

    @property
    def acquisition(self):
        return self.session.query(Candidate).filter_by(snr_id=self.snr_id, priority=0).all()

    @property
    def snr_locations(self):
        return self.session.query(SNRLocation).filter_by(snr_id=self.snr_id).all()

    @property
    def calibration_stars(self):
        return self.session.query(Candidate).filter_by(snr_id=self.snr_id, priority=3, selected=True).all()

    @property
    def all_stars(self):
        return self.session.query(Candidate).filter_by(snr_id=self.snr_id, selected=True).all()

    @misc.lazyproperty
    def mcps2cand_id(self):
        return dict([(cand.mcps.id, cand.id) for cand in self.all_stars])


    @misc.lazyproperty
    def cand2mdf_id(self):
        label = [item.label for item in self.candidates]
        mdf_id = [np.where(self.mdf['ID']==item.mcps_id)[0][0] for item in self.candidates]
        return dict(zip(label, mdf_id))

    @misc.lazyproperty
    def cand2cand_id(self):
        return dict([(item.label, item.id) for item in self.candidates])

    @misc.lazyproperty
    def mcps_id2cand_id(self):
        return dict([(item.mcps_id, item.label) for item in self.candidates])



    @misc.lazyproperty
    def slitlets_simple(self):
        slitlet_dict = {}

        for band in ('red', 'blue'):
            slitlet_dict[band] = {}
            for cw in cent_waves[band]:
                slitlet_dict[band][cw] = {}
                for raw_fits in self.science:
                    if raw_fits.header['grating'] == grating_names[band] and abs(raw_fits.header['centwave']-cw) < 1e-5:
                        fname = os.path.join(self.reduced_dir, self.snr.name, band, self.reduction_step_dict['simple'] + raw_fits.fname)
                for label in self.cand2mdf_id:
                    cur_fits = fits.open(fname)
                    slitlet_dict[band][cw][label] = cur_fits[self.cand2mdf_id[label]+2]
        return slitlet_dict

    @misc.lazyproperty
    def slitlets_wavecal(self):
        slitlet_dict = {}

        for band in ('red', 'blue'):
            slitlet_dict[band] = {}
            for cw in cent_waves[band]:
                slitlet_dict[band][cw] = {}
                for raw_fits in self.science:
                    if raw_fits.header['grating'] == grating_names[band] and abs(raw_fits.header['centwave']-cw) < 1e-5:
                        fname = os.path.join(self.reduced_dir, self.snr.name, band, self.reduction_step_dict['wavecal'] + raw_fits.fname)
                        print fname
                for label in self.cand2mdf_id:
                    cur_fits = fits.open(fname)
                    slitlet_dict[band][cw][label] = cur_fits[self.cand2mdf_id[label]+2]
        return slitlet_dict

    @misc.lazyproperty
    def slitlets_skysub(self):
        slitlet_dict = {}

        for band in ('red', 'blue'):
            slitlet_dict[band] = {}
            for cw in cent_waves[band]:
                slitlet_dict[band][cw] = {}
                for raw_fits in self.science:
                    if raw_fits.header['grating'] == grating_names[band] and abs(raw_fits.header['centwave']-cw) < 1e-5:
                        fname = os.path.join(self.reduced_dir, self.snr.name, band, self.reduction_step_dict['skysub'] + raw_fits.fname)
                        print fname
                for label in self.cand2mdf_id:
                    cur_fits = fits.open(fname)
                    slitlet_dict[band][cw][label] = cur_fits[self.cand2mdf_id[label]+2]
        return slitlet_dict


    @misc.lazyproperty
    def spectra_extracted(self):
        slitlet_dict = {}

        for band in ('red', 'blue'):
            slitlet_dict[band] = {}
            for cw in cent_waves[band]:
                slitlet_dict[band][cw] = {}
                for raw_fits in self.science:
                    if raw_fits.header['grating'] == grating_names[band] and abs(raw_fits.header['centwave']-cw) < 1e-5:
                        fname = os.path.join(self.reduced_dir, self.snr.name, band, self.reduction_step_dict['extracted'] + raw_fits.fname)
                        print fname
                for label in self.cand2mdf_id:
                    cur_fits = fits.open(fname)
                    cur_slitlet = cur_fits[self.cand2mdf_id[label]+2]
                    wave = cur_slitlet.header['crval1'] + np.arange(cur_slitlet.header['naxis1'])*cur_slitlet.header['cdelt1']
                    spec = cur_slitlet.data
                    slitlet_dict[band][cw][label] = oned.onedspec(wave, spec, mode='waveflux')

        return slitlet_dict

    def iter_stellar_spectra(self):
        snr_observations = []
        for raw_fits in self.science:
            print "Working on file %s" % raw_fits.fname
            central_wave = int(raw_fits.header['centwave'])
            band = 'red' if central_wave > 600 else 'blue'
            path =  os.path.join(self.reduced_dir, self.snr.name, band)
            fname = os.path.join(path, self.reduction_step_dict['extracted'] + raw_fits.fname)
            if not os.path.exists(fname):
                print "WARNING: %s does not exist SKIPPING" % (fname, )
                continue
            current_fits = fits.open(fname)
            current_mdf = pd.DataFrame(current_fits[1].data)
            current_mdf.set_index('ID', inplace=True)
            for fits_idx, extracted_spec in enumerate(current_fits[2:]):
                fits_idx += 2
                object_id = int(extracted_spec.header['object'])
                mdf_entry = current_mdf.ix[object_id]

                if int(mdf_entry['priority']) in (0, 1, 3):

                    #print "Found Stellar Candidate_id"
                    candidate_id = self.mcps2cand_id[object_id]
                    snr_observations.append((candidate_id, raw_fits.fname, path, central_wave, fits_idx))
                else:
                    continue

        return snr_observations














    @property
    def mask_name(self):
        return self.science[0].header['maskname']

    @property
    def mdf(self):
        if self.mdf_dir is None:
            print "mdf_dir not set"
            return None
        if self._mdf is None:
            mdf_fname = os.path.join(self.mdf_dir, self.mask_name + '.fits')
            print mdf_fname
            self._mdf = fits.getdata(mdf_fname)

        return self._mdf

class SNRNeighbour(Base):
    __tablename__ = 'snr_neighbour'
    id = Column(Integer, primary_key=True)
    mcps_id = Column(Integer, ForeignKey('mcps.id'))
    snr_id = Column(Integer, ForeignKey('snrs.id'))
    separation = Column(Float)
    mcps = relationship('MCPS', backref='separation')
    snr = relationship('SNRS', backref='snr')


    def __init__(self, mcps_id, snr_id, separation):
        self.mcps_id = mcps_id
        self.snr_id = snr_id
        self.separation = separation


class MCObject(Base):
    __tablename__ = 'mcobject'
    id = Column(Integer, ForeignKey('mcps.id'), primary_key=True)
    galex_id = Column(Integer, ForeignKey('galex.id'), default=None)
    uit_id = Column(Integer, ForeignKey('uit.id'), default=None)
    galex = relationship("Galex")
    mcps = relationship("MCPS")

    def __init__(self, galex_id, uit_id):
        self.galex_id = galex_id
        self.uit_id = uit_id


class Candidate(Base):
    __tablename__ = 'candidates'

    id = Column(Integer, primary_key=True)
    snr_id = Column(Integer, ForeignKey('snrs.id'))
    mcps_id = Column(Integer, ForeignKey('mcps.id'))
    snr_neighbour_id = Column(Integer, ForeignKey('snr_neighbour.id'))
    galex_id = Column(Integer, ForeignKey('galex.id'))
    label = Column(String)
    priority = Column(Integer)
    image_x = Column(Float)
    image_y = Column(Float)
    slit_size_x = Column(Float)
    slit_size_y = Column(Float)
    slit_tilt = Column(Float)
    slit_pos_y = Column(Float)
    selected = Column(Boolean)


    #relationship()
    galex = relationship("Galex")
    mcps = relationship("MCPS")
    snr = relationship("SNRS")
    snr_neighbour = relationship("SNRNeighbour")

    def get_ra_sex(self):
        return ephem.hours(np.deg2rad(self.ra))
    def get_dec_sex(self):
        return ephem.degrees(np.deg2rad(self.dec))

    def set_ra_sex(self, value):
        self.ra = ephem.hours(value)/ephem.degree
    def set_dec_sex(self, value):
        self.dec = ephem.degrees(value)/ephem.degree

    ra_sex = property(get_ra_sex, set_ra_sex)
    dec_sex = property(get_dec_sex, set_dec_sex)


    def __init__(self, snr_id, mcps_id, snr_neighbour_id, galex_id, label, priority, image_x=None, image_y=None, slit_size_x=1.5, slit_size_y=5, slit_tilt=0, slit_pos_y=0):
        self.snr_id = snr_id
        self.mcps_id = mcps_id
        self.snr_neighbour_id = snr_neighbour_id
        self.galex_id = galex_id
        self.label = label
        self.image_x = image_x
        self.image_y = image_y
        self.priority = priority
        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_tilt = slit_tilt
        self.slit_pos_y = slit_pos_y

    def __repr__(self):
        return "SNR Candidate MCPS #%d (%s) of SNR %s RA=%s DEC=%s V=%.2f priority=%d slit_width=%.2f slit_length=%.2f" %\
           (self.mcps_id, self.label, self.snr.name.upper(), self.mcps.ra_sex, self.mcps.dec_sex,
            self.mcps.v, self.priority, self.slit_size_x, self.slit_size_y)

class SNRLocation(Base):
    __tablename__ = 'snr_location'

    id = Column(Integer, primary_key=True)
    snr_id = Column(Integer, ForeignKey('snrs.id'))
    ra = Column(Float)
    dec = Column(Float)
    priority = Column(Integer)
    image_x = Column(Integer)
    image_y = Column(Integer)
    slit_size_x = Column(Float)
    slit_size_y = Column(Float)
    slit_tilt = Column(Float)
    slit_pos_y = Column(Float)
    selected = Column(Boolean)
    snr = relationship("SNRS")

    def get_ra_sex(self):
        return ephem.hours(np.deg2rad(self.ra))
    def get_dec_sex(self):
        return ephem.degrees(np.deg2rad(self.dec))

    def set_ra_sex(self, value):
        self.ra = ephem.hours(value)/ephem.degree
    def set_dec_sex(self, value):
        self.dec = ephem.degrees(value)/ephem.degree

    ra_sex = property(get_ra_sex, set_ra_sex)
    dec_sex = property(get_dec_sex, set_dec_sex)
    def __init__(self, snr_id, ra, dec, priority, image_x=None, image_y=None, slit_size_x=1.5, slit_size_y=5, slit_tilt=0, slit_pos_y=0):
        self.snr_id = snr_id
        self.ra = ra
        self.dec = dec
        self.image_x = image_x
        self.image_y = image_y

        self.priority = priority
        self.slit_size_x = slit_size_x
        self.slit_size_y = slit_size_y
        self.slit_tilt = slit_tilt
        self.slit_pos_y = slit_pos_y



    def __repr__(self):
        return "SNR Location #%d RA=%s DEC=%s priority=%d" % \
            (self.id, self.ra_sex, self.dec_sex, self.priority)


class MaskLabel(Base):
    __tablename__ = 'mask_label'
    id = Column(Integer, primary_key=True)
    mask_id = Column(Integer)
    snr_id = Column(Integer, ForeignKey('snrs.id'))

    snr = relationship('SNRS', backref = 'mask_label')

    def __init__(self, mask_id, snr_id):
        self.mask_id = mask_id
        self.snr_id = snr_id


class Quasar(Base):
    __tablename__ = 'quasars'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    ra = Column(Float)
    dec = Column(Float)
    v = Column(Float)
    i = Column(Float)

    in_snr_id = Column(Integer, ForeignKey('snrs.id'))
    snr_distance = Column(Float)
    mcps_id = Column(Integer, ForeignKey('mcps.id'))

    in_snr = relationship(SNRS, backref='quasars')



###### OBSERVATIONAL TABLES #########


class ObsType(Base):

    __tablename__ = 'obs_type'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    desc = Column(String)


    def __init__(self, name, desc):
        self.name = name
        self.desc = desc


class FITSFiles(object):
    @misc.lazyproperty
    def full_path(self):
        return os.path.join(self.path, self.fname)

    @misc.lazyproperty
    def data(self):
        if self._data is None:
            self._data = fits.getdata(self.full_path)
        return self._data

    @misc.lazyproperty
    def header(self):
        if self._header is None:
            self._header = fits.getheader(self.full_path)
        return self._header







class RawFITSFiles(Base, FITSFiles):
    __tablename__ = 'raw_fits_files'

    id = Column(Integer, primary_key=True)
    fname = Column(String)
    path = Column(String)
    data_label = Column(String)
    obs_type_id = Column(Integer, ForeignKey('obs_type.id'))
    obs_class = Column(String)

    object_name = Column(String)
    mask_id = Column(Integer)
    date_obs = Column(DateTime)


    snr_id = Column(Integer, ForeignKey('snrs.id'))

    snr = relationship("SNRS", backref = 'fits')


    obs_type_dict = {'object':1,
                     'flat':2,
                     'arc':3,
                     'bias':4,
                     'dark':5}
    _header = None
    _data = None

    def __init__(self, fname, path):
        self.fname = fname
        self.path = path

        header = fits.getheader(os.path.join(path, fname))

        self.data_label = header['datalab']
        self.obs_type_id = self.obs_type_dict[header['obstype'].lower()]
        self.object_name = header['object'].lower()
        self.obs_class = header['obsclass']
        self.mask_id = header['maskid']
        dateobs_iso = time.Time(header['obsepoch'], format='jyear', scale='utc').isot
        self.dateobs = parser.parse(dateobs_iso)

        self.snr_id = None


class CandidateObservation(Base):
    __tablename__ = 'candidate_observation'

    id = Column(Integer, primary_key=True)
    candidate_id = Column(Integer, ForeignKey('candidates.id'))
    snr_id = Column(Integer, ForeignKey('snrs.id'))
    fname = Column(String)
    path = Column(String)
    centralwave = Column(Integer)
    fits_idx = Column(Integer)

    candidate = relationship('Candidate', uselist=False, backref='observations')
    snr = relationship('SNRS', uselist=False, backref='observations')

    def __init__(self, candidate_id, snr_id, fname, path, centralwave, fits_idx):
        self.candidate_id = candidate_id
        self.snr_id = snr_id
        self.fname = fname
        self.path = path
        self.centralwave = centralwave
        self.fits_idx = fits_idx

    @misc.lazyproperty
    def fullpath(self):
        return os.path.join(self.path, self.fname)

    @property
    def spectrum(self):
        full_path = os.path.join(self.path, reduction_step_dict['extracted'] + self.fname)
        full_path_mask = os.path.join(self.path, reduction_step_dict['extracted'] + self.fname[:-5] + '-%d.np' % self.fits_idx)

        current_header = fits.getheader(full_path, ext=self.fits_idx)
        if os.path.exists(full_path_mask):
            mask = np.fromfile(full_path_mask,dtype=bool)
        else:
            mask = None
        wave = current_header['crval1'] + np.arange(current_header['naxis1'])* current_header['cdelt1']
        flux = fits.getdata(full_path, ext=self.fits_idx)


        return oned.onedspec(wave, flux, mask=mask, mode='waveflux')

    def __repr__(self):
        return "<Observation %s CentralWL=%d nm fits_idx = %d>" % (self.fname, self.centralwave, self.fits_idx)



class SimpleFit(Base):
    __tablename__ = 'simple_fit'

    id = Column(Integer, ForeignKey('candidate_observation.id'), primary_key=True)
    teff_bv = Column(Float)
    teff_vk = Column(Float)
    npol = Column(Integer)
    teff = Column(Float)
    logg = Column(Float)
    feh = Column(Float)

    v_rad = Column(Float)
    success = Column(Boolean)

    observation = relationship('CandidateObservation', uselist=False, backref=backref('simple_fit', uselist=False))

    @classmethod
    def from_model(cls, id, teff_bv, teff_vk, spectrum_model):
        fit = cls(teff_bv, teff_vk,
            spectrum_model.npol,
            spectrum_model.params['teff'],
            spectrum_model.params['logg'],
            spectrum_model.params['feh'],
            spectrum_model.params['v_rad'],
            spectrum_model.fit_result['success'])

        fit.id = id
        return fit

    def __init__(self, teff_bv, teff_vk, npol, teff, logg, feh, v_rad, success):
        self.teff_bv = teff_bv
        self.teff_vk = teff_vk
        self.npol = npol
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.v_rad = v_rad
        self.success = success


    def to_model(self, spectrum_model):
        spectrum_model.params['offset'] = self.offset
        spectrum_model.params['teff'] = self.teff
        spectrum_model.params['logg'] = self.logg
        spectrum_model.params['feh'] = self.feh
        spectrum_model.params['v_rad'] = self.v_rad
        spectrum_model.params['v_rot'] = self.v_rot

        return spectrum_model

class GridFit(Base):
    __tablename__ = 'grid_fit'


    id = Column(Integer, ForeignKey('candidate_observation.id'), primary_key=True)
    v_rad = Column(Float)
    teff = Column(Float)
    logg = Column(Float)
    h5fname = Column(String)

    observation = relationship('CandidateObservation', uselist=False, backref=backref('grid_fit', uselist=False))

    def __init__(self, v_rad, teff, logg, h5fname):
        self.v_rad = v_rad
        self.teff = teff
        self.logg = logg
        self.h5fname = h5fname



    @misc.lazyproperty
    def teffs(self):
        with h5py.File(self.h5fname) as h5file:
            return np.array(h5file['teffs'])

    @misc.lazyproperty
    def loggs(self):
        with h5py.File(self.h5fname) as h5file:
            return np.array(h5file['loggs'])

    @misc.lazyproperty
    def v_rads(self):
        with h5py.File(self.h5fname) as h5file:
            return np.array(h5file['v_rads'])



    @misc.lazyproperty
    def teff_logg_grid(self):
        with h5py.File(self.h5fname) as h5file:
            return np.array(h5file['teff_logg_grid'])


    @misc.lazyproperty
    def v_rad_grid(self):
        with h5py.File(self.h5fname) as h5file:
            return np.array(h5file['v_rad_grid'])

    @misc.lazyproperty
    def teff_min(self):
        logg_idx, teff_idx = np.unravel_index(np.argmin(self.teff_logg_grid), self.teff_logg_grid.shape)
        return self.teffs[teff_idx]


    @misc.lazyproperty
    def logg_min(self):
        logg_idx, teff_idx = np.unravel_index(np.argmin(self.teff_logg_grid), self.teff_logg_grid.shape)
        return self.loggs[logg_idx]


    def plot_teff_logg(self, fig):
        fig.clf()
        ax = fig.add_subplot(111)
        #grid = ax.imshow(self.teff_logg_grid.transpose(),  aspect='auto', extent=(self.teffs.min(), self.teffs.max(), self.loggs.min(), self.loggs.max()))
        #fig.colorbar(grid)
        ax.set_xlabel('Teff')
        ax.set_ylabel('log(g)')
        X, Y = np.meshgrid(self.teffs, self.loggs)
        contours = ax.contourf(X, Y, self.teff_logg_grid)
        fig.colorbar(contours)
        ax.plot([self.teff], [self.logg], 'rx', mew=3, ms=20, label='polyfit')
        ax.plot([self.teff_min], [self.logg_min], 'gx', mew=3, ms=20, label='min')
        ax.set_title('teff = %.2f teff_min = %.2f logg = %.2f logg_min = %.2f' %
                     (self.teff, self.teff_min, self.logg, self.logg_min))
        ax.legend(numpoints=1)

    def plot_vrad_chi2(self, fig):
        fig.clf()
        ax = fig.add_subplot(111)
        ax.plot(self.v_rads, self.v_rad_grid)
        ax.set_xlabel('v_rad [km/s]')
        ax.set_ylabel('Chi2')
        ax.axvline(self.v_rad, color='red', linewidth=2)
        ax.set_title('v_rad = %.2f' % self.v_rad)

    def plot_comparison(self, spectrum_model, fig):
        fig.clf()
        spectrum_model.params['teff'] = self.teff_min
        spectrum_model.params['logg'] = self.logg_min
        spectrum_model.params['feh'] = 0.0
        spectrum_model.params['v_rad'] = self.v_rad
        spectrum_model.eval_current()
        spectrum_model.plot(fig)





    def save_images(self, fig=None, spectrum_model=None):
        base_image_path = '/home/wkerzend/public_html/mcsnr/data/media/simple_spectra.d'
        image_path = os.path.join(base_image_path, self.observation.candidate.snr.name)
        image_fname = 'cand%s-%s_teff_logg_grid.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        self.plot_teff_logg(fig)
        fig.savefig(os.path.join(image_path, image_fname).lower())
        image_fname = 'cand%s-%s_v_rad_grid.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        self.plot_vrad_chi2(fig)
        fig.savefig(os.path.join(image_path, image_fname).lower())
        image_fname = 'cand%s-%s_grid_comparison.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        self.plot_comparison(spectrum_model, fig)
        fig.savefig(os.path.join(image_path, image_fname).lower())

    def make_wiki_table(self, image_width=600):

        image_template = '{{simple_spectra.d:%(snr_name)s:%(fname)s?%(width)s}}'
        wiki_str = "^ Teff-Logg-Grid      ^ V rad grid       ^ Comparison spectra       ^\n"

        image_fname = 'cand%s-%s_teff_logg_grid.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        current_wiki_image = image_template % dict(snr_name=self.observation.snr.name, fname=image_fname, width=image_width)
        wiki_str += "|%s" % current_wiki_image

        image_fname = 'cand%s-%s_v_rad_grid.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        current_wiki_image = image_template % dict(snr_name=self.observation.snr.name, fname=image_fname, width=image_width)
        wiki_str += "|%s" % current_wiki_image

        image_fname = 'cand%s-%s_grid_comparison.png' % (self.observation.candidate.label, self.observation.fname[:-5])
        current_wiki_image = image_template % dict(snr_name=self.observation.snr.name, fname=image_fname, width=image_width)
        wiki_str += "|%s|\n" % current_wiki_image

        return wiki_str








def make_model(self, spectrum_model, teff, logg, feh=0.0):
            model = spectrum_model()



class MartenFit(Base):
    __tablename__ = 'marten_fit'

    id = Column(Integer, ForeignKey('candidate_observation.id'), primary_key=True)
    teff = Column(Float)
    p0 = Column(Float)
    p1 = Column(Float)
    p2 = Column(Float)
    p3 = Column(Float)
    p4 = Column(Float)
    p5 = Column(Float)
    v_rad = Column(Float)
    v_rad_error = Column(Float)

    observation = relationship('CandidateObservation', uselist=False, backref=backref('marten_fit', uselist=False))


    def __init__(self, teff, p0, p1, p2, p3, p4, p5, v_rad, v_rad_error):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.p5 = p5
        self.teff = teff
        self.v_rad = v_rad
        self.v_rad_error = v_rad_error

    def to_model(self, spectrum_model):
        spectrum_model.params['p0'] = self.p0
        spectrum_model.params['p1'] = self.p1
        spectrum_model.params['p2'] = self.p2
        spectrum_model.params['p3'] = self.p3
        spectrum_model.params['p4'] = self.p4
        spectrum_model.params['offset'] = self.offset
        spectrum_model.params['teff'] = self.teff
        spectrum_model.params['logg'] = self.logg
        spectrum_model.params['feh'] = self.feh
        spectrum_model.params['v_rad'] = self.v_rad
        spectrum_model.params['v_rot'] = self.v_rot

        return spectrum_model






def init_obs_type(session):

    obstypes = []

    obs_type_dict = {'object':1,
                     'flat':2,
                     'arc':3,
                     'bias':4,
                     'dark':5}
    for key, item in obs_type_dict.items():
        new_obs_type = ObsType(key, None)
        new_obs_type.id = item
        obstypes.append(new_obs_type)

    session.add_all(obstypes)
    session.commit()



def validate_mag(mag, mag_err, lower_limit=2, upper_limit=30):
    if lower_limit <= mag <= upper_limit:
        return mag, mag_err
    else:
        return None, None

def initialize_support_tables():
    pass

def initialize_session(engine):
    Base.metadata.bind=engine
    Session = sessionmaker(bind=engine)
    session = Session()
    return session

def update_mcps(session, min_mag=2, max_mag=30):
    mcps = MCPS.__table__
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.u,min_mag, max_mag))).values(u=None, u_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.b,min_mag, max_mag))).values(b=None, b_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.v,min_mag, max_mag))).values(v=None, v_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.i,min_mag, max_mag))).values(i=None, i_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.j,min_mag, max_mag))).values(j=None, j_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.h,min_mag, max_mag))).values(h=None, h_err=None))
    session.execute(update(mcps).where(sql.not_(sql.between(mcps.c.k,min_mag, max_mag))).values(k=None, k_err=None))
    session.commit()



def ingest_uit(session, engine):
    UIT.__table__.drop(engine, checkfirst=True)
    UIT.__table__.create(engine)
    smc_raw = np.recfromtxt('uit/smc.dat',
        usecols= (4,5,6,7,8,9, 10, 11),
        names=('rah', 'ram', 'ras', 'decd', 'decm', 'decs', 'u', 'u_err'))
    ras= 15*(smc_raw['rah'] + smc_raw['ram'] / 60. + smc_raw['ras'] / 3600.)
    decs = smc_raw['decd'] + smc_raw['decm'] / 60. + smc_raw['decs'] / 3600.
    for ra, dec, m162, m162_err in zip(ras, decs, smc_raw['u'], smc_raw['u_err']):
        session.add(UIT(ra, dec, m162, m162_err, None, None, None, None, 2))


    lmc_raw = np.recfromtxt('uit/lmc.dat',
        usecols=(1, 2, 3, 4, 5, 6, 9, 10, 12, 13),
        names=('rah', 'ram', 'ras', 'decd', 'decm', 'decs', 'b1', 'b1_err', 'b5', 'b5_err'))
    ras= 15*(lmc_raw['rah'] + lmc_raw['ram'] / 60. + lmc_raw['ras'] / 3600.)
    decs = lmc_raw['decd'] + lmc_raw['decm'] / 60. + lmc_raw['decs'] / 3600.
    for ra, dec, b1, b1_err, b5, b5_err in zip(ras, decs, lmc_raw['b1'], lmc_raw['b1_err'], lmc_raw['b5'], lmc_raw['b5_err']):
        session.add(UIT(ra, dec, None, None, b1, b1_err, b5, b5_err, 1))


    session.commit()

def ingest_galex(session, engine):

    GalexBand.__table__.drop(engine, checkfirst=True)
    GalexBand.__table__.create(engine)
    session.add(band1)
    session.add(band2)
    session.add(band3)
    session.commit()
    Galex.__table__.drop(engine, checkfirst=True)
    Galex.__table__.create(engine)
    smc_raw = np.recfromcsv('galex_db/smc_wkerzendorf.csv', names=True)
    for line in smc_raw:
        session.add(Galex(line['ra'], line['dec'],
            line['nuv_mag'], line['nuv_magerr'],
            line['fuv_mag'], line['fuv_magerr'],
            line['fov_radius'], int(line['band']), 2))

    lmc_raw = np.recfromcsv('galex_db/lmc_wkerzendorf.csv', names=True)
    for line in lmc_raw:
        session.add(Galex(line['ra'], line['dec'],
            line['nuv_mag'], line['nuv_magerr'],
            line['fuv_mag'], line['fuv_magerr'],
            line['fov_radius'], int(line['band']), 1))
    session.commit()


def get_session():
    pass

def get_conn():
    pass
