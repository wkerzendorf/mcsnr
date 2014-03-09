from geminiutil.gmos import GMOSMOSProject
from geminiutil.gmos import gmos_alchemy as ga


import socket
host = socket.gethostname()
if host == 'swan':
    work_dir = '/work/mhvk/mcsnr/'
elif host == 'moria':
    work_dir = '/media/data1/mcsnr/'
else:
    raise ValueError("Unknown host='{0}'".format(host))

dbname = 'sqlite:///gmos.db3'

proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)
