C:
cd C:\Program Files\FlightGear

SET FG_ROOT=C:\Program Files\FlightGear\\data
.\\bin\\win64\\fgfs --aircraft=f16 --fdm=network,localhost,5501,5502,5503 --fog-fastest --disable-clouds --start-date-lat=2004:06:01:09:00:00 --disable-sound --in-air --enable-freeze --airport=EGTC --runway=1 --altitude=14435.7 --heading=0 --offset-distance=4.72 --offset-azimuth=0
