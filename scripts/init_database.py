from geminiutil.gmos import GMOSMOSProject
import os
import socket
host = socket.gethostname()
if host == 'swan':
    work_dir = '.'
    dbname = 'sqlite:///gmos.db3'
elif host == 'moria':
    work_dir = '/media/data1/mcsnr/'
    dbname = 'sqlite:///databases/gmos_working_copy.db3'
else:
    raise ValueError("Unknown host='{0}'".format(host))


os.system('rm {0}'.format(db_fname))

gmos_working_copy.db3

proj = GMOSMOSProject(dbname, work_dir=work_dir, echo=False)
proj.initialize_database()
proj.add_directory('gmos_data/raw')
proj.add_directory('gmos_data/mdf_dir')
proj.link_database()
