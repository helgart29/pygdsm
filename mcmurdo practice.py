import numpy as np
import pylab as plt
import healpy as hp
import os

from pygdsm import LFSMObserver, LowFrequencySkyModel, GlobalSkyModel, GSMObserver
from datetime import datetime, timezone
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

(latitude, longitude, elevation) = ('-77.846', '166.676', 183)

ov = GSMObserver()
ov.lon = longitude
ov.lat = latitude
ov.elev = elevation
ov.date = datetime(2025, 12, 1, 0, 0)
ov.generate(150)
frequency=150
d = ov.view_observed_gsm(logged=True, show=False)

plt.savefig('mcmurdo_mollweide_updated.png')


#attempt at writing mollweide
def observed_gsm(self):
        """ Return the GSM (Mollweide), with below-horizon area masked. """
        sky = self.observed_sky

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi / 2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        derotate = hp.Rotator(rot=[ra_deg, dec_deg])
        g0, g1 = derotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        coordrotate = hp.Rotator(coord=['C', 'G'], inv=True)
        g0, g1 = coordrotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]
        return sky

def updated_observed_gsm():
    (latitude, longitude, elevation) = ('-77.846', '166.676', 183)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation

    ov.date=datetime(2025,12,1,0,0)
    ov.generate(150)
    
    obs = []

    sky = ov.view_observed_gsm(logged=True, show=False)

    #hp.mollview(sky, coord=['G','C'], title='Mollweide View from McMurdo Station at 150 MHz')
    #hp.graticule(coord=['C'])
    hp.projview(sky, coord = ['G','C'], unit="Flux Density[Jy]",fontsize=dict(title='large',
                  xlabel='medium', ylabel='medium', 
                  cbar_label='small'),
                  graticule=True, graticule_labels=True, projection_type='mollweide',
                  title ='Mollweide View from McMurdo Station at 150 MHz', xlabel='Right Ascension', ylabel='Declination')
    plt.savefig('testing updated observer.png')



def galaxy_mcmurdo_check():
    (latitude, longitude, elevation) = ('-77.846', '166.676', 183)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation

    ov.date=datetime(2025,12,1,0,0)
    ov.generate(150)
    
    obs = []

    sky = ov.view_observed_gsm(logged=True, show=False)

    #hp.mollview(sky, coord=['G','C'], title='Mollweide View from McMurdo Station at 150 MHz')
    #hp.graticule(coord=['C'])
    hp.projview(sky, coord = ['G','C'], unit="Flux Density[Jy]",fontsize=dict(title='large',
                  xlabel='medium', ylabel='medium', 
                  cbar_label='small'),
                  graticule=True, graticule_labels=True, projection_type='mollweide',
                  title ='Mollweide View from McMurdo Station at 150 MHz', xlabel='Right Ascension', ylabel='Declination')
    
    #ra_point = 266.417
    #dec_point = -29.00781

    #theta_point = np.deg2rad(90 - dec_point)
    #phi_point = np.deg2rad(ra_point)

    hp.projtext(266.417, -29.00781, lonlat=True, coord='C')
    
    plt.savefig('testing galaxy check.png')

galaxy_mcmurdo_check()