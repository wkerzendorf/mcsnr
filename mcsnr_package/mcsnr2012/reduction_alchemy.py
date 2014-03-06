__author__ = 'wkerzend'





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
