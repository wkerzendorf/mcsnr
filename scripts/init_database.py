from geminiutil.gmos import GMOSMOSProject


dbname = 'sqlite:///databases/gmos.db3'

proj = GMOSMOSProject(dbname, work_dir='/media/data1/mcsnr/', echo=False)

proj.initialize_database()
proj.add_directory('gmos_data/raw')
proj.add_directory('gmos_data/mdf_dir')
proj.link_database()
