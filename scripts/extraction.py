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
