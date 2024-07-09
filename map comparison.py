import numpy as np
import pylab as plt
import healpy as hp
import os
import inspect

from pygdsm import GSMObserver, GSMObserver16, LFSMObserver, GlobalSkyModel
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from datetime import datetime
import astropy.units as u


# (latitude, longitude, elevation) = ('-77.853','167.213',-0.01097561)
# ov = GSMObserver()
# ov.lon = longitude
# ov.lat = latitude
# ov.elevation = elevation
# ov.date = datetime(2016,12,2,13,11,0)
# ov.generate(50)
# sky = ov.view_observed_gsm(logged=True, show=False)
# hp.projview(sky, coord = ['G','C'], unit="log\N{SUBSCRIPT TWO}(T) [K]",fontsize=dict(title='large',
#                   xlabel='medium', ylabel='medium', 
#                   cbar_label='small'),
#                   graticule=True, graticule_labels=True, projection_type='mollweide',
#                   title = f'GSM08 McMurdo Station at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
# plt.savefig('GSM08 McMurdo')


(latitude, longitude, elevation) = ('-77.853','167.213',-0.01097561)
ov = GSMObserver16()
ov.lon = longitude
ov.lat = latitude
ov.elevation = elevation
ov.date = datetime(2016,12,2,13,11,0)
ov.generate(50)
sky = ov.view_observed_gsm()
hp.projview(sky, coord = ['G','C'], unit="Temperature [K]",fontsize=dict(title='large',
                  xlabel='medium', ylabel='medium', 
                  cbar_label='small'),min = 500, max = 60000,
                  graticule=True, graticule_labels=True, projection_type='mollweide',
                  title = f'GSM16 McMurdo Station at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
# plt.savefig('GSM16 McMurdo')

# (latitude, longitude, elevation) = ('-77.853','167.213',-0.01097561)
# ov = LFSMObserver()
# ov.lon = longitude
# ov.lat = latitude
# ov.elevation = elevation
# ov.date = datetime(2016,12,2,13,11,0)
# ov.generate(50)
# sky = ov.view_observed_gsm(logged=True)
# hp.projview(sky, coord = ['G','C'], unit="log\N{SUBSCRIPT TWO}(T) [K]",fontsize=dict(title='large',
#                   xlabel='medium', ylabel='medium', 
#                   cbar_label='small'),
#                   graticule=True, graticule_labels=True, projection_type='mollweide',
#                   title = f'LFSM McMurdo Station at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
# plt.savefig('LFSM McMurdo')


def test_get_sky_temperature():
    gc = SkyCoord(266.41683, -29.00781, unit='deg', frame='icrs')
    frequency = 100
    g = GlobalSkyModel()
    T = g.get_sky_temperature(gc, frequency, include_cmb=True)
    
print(test_get_sky_temperature())

