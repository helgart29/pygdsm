import numpy as np
import h5py
import healpy as hp
from astropy.io import fits
from .plot_utils import show_plt

def is_fits(filepath):
    """
    Check if file is a FITS file
    Returns True of False

    Parameters
    ----------
    filepath: str
        Path to file
    """
    FITS_SIGNATURE = (b"\x53\x49\x4d\x50\x4c\x45\x20\x20\x3d\x20\x20\x20\x20\x20"
                      b"\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20\x20"
                      b"\x20\x54")
    with open(str(filepath),'rb') as f:
        try:
            return f.read(30) == FITS_SIGNATURE
        except FileNotFoundError as e:
            print(e)
            return False


class BaseSkyModel(object):
    """ Global sky model (GSM) class for generating sky models.
    """
    def __init__(self, name, filepath, freq_unit, data_unit, basemap):
        """ Initialise basic sky model class

        Parameters
        ----------
        name (str):      Name of GSM
        filepath (str):    Path to HDF5 data / FITS data (healpix) to load
        freq_unit (str): Frequency unit (MHz / GHz / Hz)
        data_unit (str): Unit for pixel scale (e.g. K)
        basemap (str):   Map used as a basis for spatial structure in PCA fit.

        Notes
        -----
        Any GSM needs to supply a generate() function
        """
        self.name = name
        if h5py.is_hdf5(filepath):
            self.h5 = h5py.File(filepath, "r")
        elif is_fits(filepath):
            self.fits = fits.open(filepath, "readonly")
        else:
            raise RuntimeError(f"Cannot read HDF5/FITS file {filepath}")
        self.basemap = basemap
        self.freq_unit = freq_unit
        self.data_unit = data_unit

        self.generated_map_data = None
        self.generated_map_freqs = None

    def generate(self, freqs):
        raise NotImplementedError

    def view(self, idx=0, logged=False, show=False):
        """ View generated map using healpy's mollweide projection.

        Parameters
        ----------
        idx: int
            index of map to view. Only required if you generated maps at
            multiple frequencies.
        logged: bool
            Take the log of the data before plotting. Defaults to False.

        """

        if self.generated_map_data is None:
            raise RuntimeError("No GSM map has been generated yet. Run generate() first.")

        if self.generated_map_data.ndim == 2:
            gmap = self.generated_map_data[idx]
            freq = self.generated_map_freqs[idx]
        else:
            gmap = self.generated_map_data
            freq = self.generated_map_freqs

        if logged:
            gmap = np.log2(gmap)

        hp.mollview(gmap, coord='G', title='%s %s, %s' % (self.name, str(freq), self.basemap))

        if show:
            show_plt()

    def write_fits(self, filename):
        """ Write out map data as FITS file.

        Parameters
        ----------
        filename: str
            file name for output FITS file
        """
        hp.write_map(filename, self.generated_map_data, column_units=self.data_unit)