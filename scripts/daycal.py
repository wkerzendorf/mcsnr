execfile("open_db.py")

daycalfits = proj.session.query(ga.GMOSLongSlitArc).all()

for daycal in daycalfits[0:1]:
    try:
        daycal.raw.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    daycal.longslit_calibrate(doplot=True)
    # daycal.longslit_calibrate_to_database()

# stores results in daycal.wave_cal.full_path
