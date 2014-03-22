from geminiutil.gmos import GMOSMOSProject

import socket
host = socket.gethostname()
if host == 'swan':
    work_dir = '.'
    dbname = 'sqlite:///gmos.db3'
elif host == 'moria':
    work_dir = '/media/data1/mcsnr/'
    dbname = 'sqlite:///databases/gmos.db3'
else:
    raise ValueError("Unknown host='{0}'".format(host))


proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)
proj.initialize_database()
proj.add_directory('gmos_data/raw')
proj.add_directory('gmos_data/mdf_dir')
proj.link_database()
