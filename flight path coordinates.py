import numpy as np
import pylab as plt
import healpy as hp
import os
import itertools
#import cv2
from PIL import Image
from moviepy.editor import *

from pygdsm import LFSMObserver, LowFrequencySkyModel, GlobalSkyModel, GSMObserver
from datetime import datetime, timezone
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

import pandas as pd
import matplotlib.animation as animation

#Getting lat and lon coordinates
data = pd.read_csv("675N (everything).csv")
print(data.columns)

df = pd.DataFrame(data)

#converting data into arrays
lat = df['Latitude (degrees)'].to_numpy()
lon = df['Longitude (degrees)'].to_numpy()
times = df['Time (seconds since launch)'].to_numpy()
elevs = df['Altitude (km)'].to_numpy()*1000

# def testing_nominal_skymap():

#     ov = GSMObserver()
    
#     obs = []

#     lats = list(range(lat))
#     lats



#     for i in lat,lon,elevs:
#         (latitude, longitude, elevation) = (lat[:i*1000],lon[:i*1000],elevs[:i*1000])

#         ov.lon = longitude
#         ov.lat = latitude
#         ov.elev = elevation


#         for i in range(times):
#             ov.date = datetime(2025,12,1,0,0,i)
#             ov.generate(50)
#             print(i)

#             sky = ov.view_observed_gsm(logged=True, show=False)

    #hp.mollview(sky, coord=['G','C'], title='Mollweide View from McMurdo Station at 150 MHz')
    #hp.graticule(coord=['C'])
    # hp.projview(sky, coord = ['G','C'], unit="Flux Density[Jy]",fontsize=dict(title='large',
    #               xlabel='medium', ylabel='medium', 
    #               cbar_label='small'),
    #               graticule=True, graticule_labels=True, projection_type='mollweide',
    #               title = f'Mollweide View from McMurdo Station at {f} MHz', xlabel='Right Ascension', ylabel='Declination')
    # plt.savefig(f'Mollweide View from McMurdo Station at {f} MHz.png')

lats = lat[::500]
lons = lon[::500]
time = times[::500]
height = elevs[::500]
label = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,47,38,39,40,41,42]
print(lats)
# print(lons)
# print(height)
# for longitude in lons:
#     ov = GSMObserver()
#     ov.lon = longitude
#     print (longitude)

# for latitude in lats:
#     ov.lat = latitude
#     print(latitude)

# for elev in height:
#     ov.elev = elev
#     print(elev)

# for (latitude, longitude, elevation, lab) in zip (lats, lons, height, label):
#     # (latitude, longitude, elevation) = ('lats', 'lons', height)
#     ov = GSMObserver()
#     ov.lon = np.array2string(longitude)
#     ov.lat = np.array2string(latitude)
#     ov.elev =  elevation

#     ov.date = datetime(2025,12,1,0,0)
#     ov.generate(50)
           
#     sky = ov.view_observed_gsm(logged=True, show=False)

#     hp.projview(sky, coord = ['G','C'], unit="Flux Density[Jy]",fontsize=dict(title='large', 
#                         xlabel='medium', ylabel='medium', cbar_label='small'),
#                         graticule=True, graticule_labels=True, projection_type='mollweide',
#                         title = f'Mollweide View from ANITA IV at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
#     plt.savefig(f'plot {lab}.png')
    


def updated_observed_gsm():
    (latitude, longitude, elevation) = ('-77.853', '167.213', -10.9756100)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation

    ov.date=datetime(2025,12,1,0,0)
    ov.generate()
    
    obs = []
    sky = ov.view_observed_gsm(logged=True, show=False)

    #hp.mollview(sky, coord=['G','C'], title='Mollweide View from McMurdo Station at 150 MHz')
    #hp.graticule(coord=['C'])
    hp.projview(sky, coord = ['G','C'], unit="Flux Density[Jy]",fontsize=dict(title='large',
                  xlabel='medium', ylabel='medium', 
                  cbar_label='small'),
                  graticule=True, graticule_labels=True, projection_type='mollweide',
                  title = f'Mollweide View from McMurdo Station at 50 MHz', xlabel='Right Ascension', ylabel='Declination')
    plt.savefig(f'plot 1.png')

clip = ImageSequenceClip(['plot 0.png', 'plot 1.png','plot 2.png','plot 3.png','plot 4.png','plot 5.png','plot 6.png','plot 7.png','plot 8.png',
                          'plot 9.png','plot 10.png','plot 11.png','plot 12.png','plot 13.png','plot 14.png','plot 15.png',
                           'plot 16.png','plot 17.png','plot 18.png','plot 19.png','plot 20.png', 'plot 21.png','plot 22.png', 'plot 23.png', 
                           'plot 24.png', 'plot 25.png', 'plot 26.png', 'plot 27.png', 'plot 28.png', 'plot 29.png', 'plot 30.png', 'plot 31.png', 
                           'plot 32.png', 'plot 33.png', 'plot 34.png', 'plot 35.png', 'plot 36.png', 'plot 47.png', 'plot 38.png', 'plot 39.png', 
                           'plot 40.png', 'plot 41.png'], fps=4)

clip.write_videofile("Galactic Emission for ANITA IV.mp4", codec="libx264")