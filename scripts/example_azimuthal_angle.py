from pysolar.solar import *
import datetime

date =datetime.datetime(2020, 1, 31, 15, 13, 1, 0, tzinfo=datetime.timezone.utc)
azim = get_azimuth(13.16, -59.2, date)

print(azim)
