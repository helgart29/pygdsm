import numpy as np
import pylab as plt
import healpy as hp
import os
import inspect
import pygdsm

from pygdsm import base_observer, GlobalSkyModel
from base_observer import BaseObserver
from gsm08 import GlobalSkyModel, GSMObserver
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from datetime import datetime
import astropy.units as u
from base_skymodel import BaseSkyModel


def updated_observed_gsm(f):
    (latitude, longitude, elevation) = ('-77.853','167.213',-10.97561)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation

    ov.date=datetime(2016,12,2,13,11)
    ov.generate(f)
    
    # obs = []
    # sky = ov.view_observed_gsm(logged=True, show=False)

  
    # #hp.mollview(sky, coord=['G','C'], title='Mollweide View from McMurdo Station at 150 MHz')
    # #hp.graticule(coord=['C'])
    # hp.projview(sky, coord = ['G','C'], unit="log\N{SUBSCRIPT TWO}(T) [K]",fontsize=dict(title='large',
    #               xlabel='medium', ylabel='medium', 
    #               cbar_label='small'),
    #               graticule=True, graticule_labels=True, projection_type='mollweide',
    #               title = f'Mollweide View from McMurdo Station at {f} MHz', xlabel='Right Ascension', ylabel='Declination')
    # plt.savefig(f'Unmasked McMurdo Station at {f} MHz.png')

updated_observed_gsm(50)


