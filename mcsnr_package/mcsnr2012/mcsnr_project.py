from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mcsnr2012 import mcsnr_alchemy as malchemy
from mcsnr2012.mcsnr_alchemy_mixins import GeminiUtilDBMixin
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import cm

from mpld3 import plugins
from mcsnr_mpld3 import css, LinkedView

from geminiutil.gmos.alchemy.mos import MOSPointSource, MOSSpectrum

from mcsnr2012.spectral_fitting import get_spectral_fit

class MCSNRProject(object):

    quick_search_classes = [malchemy.Candidate]

    quick_search_classes = [malchemy.Candidate]
    def __init__(self, database_string, echo=False):
        self.metadata = malchemy.Base.metadata
        self.engine = create_engine(database_string, echo=echo)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.conn = self.session.bind.connect()

        MOSSpectrum.spectral_grid = None



    @property
    def quick_search_table_names(self):
        return [class_model.__tablename__ for class_model in self.quick_search_classes]

    @property
    def observed_snrs(self):
        return self.session.query(malchemy.SNRGeminiTarget).filter_by(observed=True)

    def set_geminiutil_session(self, geminiutil_session):
        GeminiUtilDBMixin.geminiutil_session = geminiutil_session

    def plot_snr_cmd(self, snr, ax, color='bv', mag='v', add_candidates=True):
        cmd_search = self.session.query(malchemy.MCPS).\
            join(malchemy.SNRNeighbour).join(malchemy.SNRS).\
            filter(malchemy.SNRS.id==snr.snr.id)

        cmd_data = [(item.__getattribute__(color[0]) -
                     item.__getattribute__(color[1]),
                     item.__getattribute__(mag)) for item in cmd_search]
        data = np.array(cmd_data)
        H, xedges, yedges = np.histogram2d(data[:,1], data[:,0], bins=100,
                                           range=[[16,22],[-1,2]])
        extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
        ax.imshow(H, extent=extent, interpolation='nearest',
                  cmap=cm.gray_r,
                  aspect='auto')

        ax.set_xlabel('Color {0} - {1}'.format(*list(color)))
        ax.set_ylabel('mag {0}'.format(mag))

        if add_candidates:
            candidate_coord = [(item.mcps.__getattribute__(color[0]) -
                     item.mcps.__getattribute__(color[1]),
                     item.mcps.__getattribute__(mag)) for item in snr.candidates]

            candidate_coord = np.array(candidate_coord)
            cand_plot = ax.plot(candidate_coord[:,0], candidate_coord[:,1], 'bo')
            cand_labels = [str(item) for item in snr.candidates]
            tooltip = plugins.PointHTMLTooltip(cand_plot[0], cand_labels,
                                   voffset=10, hoffset=10)

            plugins.connect(ax.figure, tooltip)



        return cand_plot

    def plot_interactive_cmd_spectrum(self, snr, fig=None, color='bv', mag='v'):
        ax = []
        if fig is None:
            fig = plt.gcf()

        ax.append(fig.add_axes([0.05, 0.05, 0.2, 0.8]))
        ax.append(fig.add_axes([0.3, 0.1, 0.65, 0.8]))

        ax = ax[::-1]
        points = self.plot_snr_cmd(snr, ax[1], color=color, mag=mag)
        spectrum_data = []
        for cand in snr.candidates:
            spectrum_data.append([cand.mos_spectra[0].wavelength.value,
                                  cand.mos_spectra[0].flux])

        spectrum_data = np.array(spectrum_data)

        spectrum_data_interactive = spectrum_data.transpose(0, 2, 1).tolist()
        lines = ax[0].plot(cand.mos_spectra[0].wavelength.value, cand.mos_spectra[0].wavelength.value * 0, '-w', lw=3, alpha=0.5)
        ax[0].set_ylim(0, 2000)
        plugins.connect(fig, LinkedView(points[0], lines[0], spectrum_data_interactive))


        return ax

