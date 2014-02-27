from geminiutil.gmos import GMOSMOSProject

import socket
host = socket.gethostname()
if host == 'swan':
    work_dir = '/work/mhvk/mcsnr/'
elif host == 'moria':
    work_dir = '/media/data1/mcsnr/'
else:
    raise ValueError("Unknown host='{0}'".format(host))

dbname = 'sqlite:///databases/gmos.db3'

proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)
proj.initialize_database()
proj.add_directory('gmos_data/raw')
proj.add_directory('gmos_data/mdf_dir')
proj.link_database()
