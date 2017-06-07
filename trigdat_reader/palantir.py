import numpy as np
import gbmgeometry
import matplotlib.pyplot as plt
from matplotlib.projections.geo import GeoAxes
import healpy as hp

import astropy.units as u


class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""
    extrapi=False

    def __init__(self,*a,**aa):
        if 'extrapi' in aa:
            self.extrapi=aa['extrapi']
        del aa['extrapi']
        super(ThetaFormatterShiftPi,self).__init__(*a,**aa)

    def __call__(self, x, pos=None):
        if self.extrapi:
            x+=np.pi
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)




def check_power_of_two(num):
    """Check that the given number is a power of 2"""

    return num != 0 and ((num & (num - 1)) == 0)


def pix_to_sky(idx, nside):
    """Convert the pixels corresponding to the input indexes to sky coordinates (RA, Dec)"""

    theta, phi = hp.pix2ang(nside, idx)

    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)

    return ra, dec


def sky_to_pix(ra, dec, nside):
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    ipix = hp.ang2pix(nside, theta, phi)
    return ipix


_det_kwargs = {'colorbar': True,
               'cmap':'viridis',
               'show_earth':True,
               'width':18,
               'fov': 60
               }


class Palantir(object):
    def __init__(self, result, nside=64, trigdat=None, time=0.):


        """

        View teh 3ML BALROG results

        :param result: 3ML BALROG result
        :param nside: power of 2
        :param trigdat: optional trigdat files
        :param time: time of observation
        """
        self._nside = nside

        if not check_power_of_two(self._nside):
            raise RuntimeError("nside must be a power of 2.")

        self._map = np.arange(hp.nside2npix(nside), dtype=float)
        self._map[:] = 0.

        ra = result.samples[0, :]
        dec = result.samples[1, :]

        for r, d in zip(ra, dec):
            self._map[sky_to_pix(r, d, nside)] += 1

        # self._map /= float(len(result.samples))
        self._map /= float(self._map.sum())

        self._trigdat = trigdat

        if trigdat is not None:

            self._position_interpolator = gbmgeometry.PositionInterpolator(trigdat=trigdat)

            self._gbm = gbmgeometry.GBM(self._position_interpolator.quaternion(time),
                                        self._position_interpolator.sc_pos(time) * u.km)

        else:

            self._position_interpolator = None

    def skymap(self, *dets, **kwargs):

        """

        :param dets: optional detectors
        :param cmap: colormap for location
        :param show_earth: show the earth points
        :param width: no idea
        :param fov: FOV of GBM detectors
        :return: figure
        """
        vmin = 0.
        vmax = self._map.max()

        xsize = 2000
        ysize = xsize / 2.

        for kw, val in kwargs.iteritems():

            _det_kwargs[kw] = val


        theta = np.linspace(np.pi, 0, ysize)


        phi = np.linspace(-np.pi, np.pi, xsize)
        longitude = np.radians(np.linspace(-180, 180, xsize))
        latitude = np.radians(np.linspace(-90, 90, ysize))

        # project the map to a rectangular matrix xsize x ysize
        PHI, THETA = np.meshgrid(phi, theta)
        grid_pix = hp.ang2pix(self._nside, THETA, PHI)

        grid_map = self._map[grid_pix]

        # for width in [ 8.8]:
        # for width in [18., 12., 8.8]:

        fig, ax = plt.subplots(subplot_kw=dict(projection='mollweide'))

        # rasterized makes the map bitmap while the labels remain vectorial
        # flip longitude to the astro convention
        image = ax.pcolormesh(longitude[::-1],
                              latitude,
                              grid_map,
                              vmin=vmin,
                              vmax=vmax,
                              rasterized=True,
                              cmap=_det_kwargs['cmap'], zorder=-32)



        # graticule
        lon_g = ax.set_longitude_grid(60)
        ax.xaxis.set_major_formatter(ThetaFormatterShiftPi(60, extrapi=False))

        width = _det_kwargs['width']

        if width < 10:
            lat_g = ax.set_latitude_grid(45)
            lon_g = ax.set_longitude_grid_ends(90)

        for g in ax.xaxis.get_gridlines() + ax.yaxis.get_gridlines():
            g.set_linestyle("dotted")
            g.set_color("black")
            g.set_alpha(0.5)

        for det in dets:
            ra, dec = self._get_detector_map(det, _det_kwargs['fov'])

            # pix = np.where(detector_map == 1)[0]
            # ra, dec = pix_to_sky(pix,self._nside)


            x = np.radians(ra - 180)
            y = np.radians(dec)

            idx = y.argsort()

            ax.plot(x[idx], y[idx], '.', markersize=4)

        if _det_kwargs['show_earth']:
            earth_points = self._gbm.get_earth_points(False)

            ax.plot(np.radians(earth_points.ra.value - 180),
                    earth_points.dec.radian,
                    '.',
                    color='k',
                    alpha=.5,
                    zorder=-31)

        unit = r'$p$'

        # colorbar
        if _det_kwargs['colorbar']:

            cb = fig.colorbar(image,
                              orientation='horizontal',
                              shrink=.6, pad=0.05,
                              ticks=[vmin, vmax])

            cb.ax.xaxis.set_label_text(unit)
            cb.ax.xaxis.labelpad = -8
            # workaround for issue with viewers, see colorbar docstring
            cb.solids.set_edgecolor("face")

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10, color='k')
        ax.grid(True)

        return fig

    def _get_detector_map(self, det,fov):

        det_map = np.arange(hp.nside2npix(self._nside), dtype=float)
        det_map[:] = 0.

        cir = self._gbm.detectors[det].get_fov(fov)
        ra = cir[0]
        dec = cir[1]

        return ra, dec














