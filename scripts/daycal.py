execfile("open_db.py")

daycalfits = proj.session.query(ga.GMOSLongSlitArc).all()
daycalfits = [fil for fil in daycalfits
              if fil.raw.mask.name.startswith('0.5')]

for daycal in daycalfits:
    try:
        daycal.raw.prepare_to_database()
    except ga.GMOSDatabaseDuplicate:
        pass

    # arctab, linesall, shift, fake = daycal.longslit_calibrate(doplot=True)
    daycal.longslit_calibrate_to_database(force=True)

# stores results in daycal.wave_cal.full_path
