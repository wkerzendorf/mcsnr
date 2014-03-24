execfile("open_db.py")

for sci_set in proj.science_sets:
    try:
        sci_set.science.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    try:
        sci_set.calculate_slice_geometries_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    try:
        sci_set.mask_arc.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    sci_set.extract_spectra_and_calibrate(priorities=[1, 2])
