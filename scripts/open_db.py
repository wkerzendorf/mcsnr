from geminiutil.gmos import GMOSMOSProject
from geminiutil.gmos import gmos_alchemy as ga


import socket
host = socket.gethostname()
if host == 'swan':
    work_dir = '/work/mhvk/mcsnr/'
    dbname = 'sqlite:///gmos.db3'

elif host == 'moria':
    work_dir = '/media/data1/mcsnr/'
    dbname = 'sqlite:///databases/gmos_working_copy.db3'

else:
    raise ValueError("Unknown host='{0}'".format(host))


gmos_working_copy.db3

proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)
