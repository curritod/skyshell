import numpy as np
from cartopy import crs

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import GeocentricTrueEcliptic

import matplotlib.pyplot as plt
plt.style.use('dark_background')

from .catalogue_helpers import split_on_wrap
from .catalogue_helpers import load_hipparcos_positions
from .catalogue_helpers import load_constellation_polylines

from .sun_moon_position import get_sun_position
from .sun_moon_position import get_moon_position

#Default observatory
from astropy.coordinates import EarthLocation
PARANAL = EarthLocation(lat=-24.6272*u.deg, lon=-70.4042*u.deg, height=2635*u.m)

#Load package resources
from importlib.resources import files

hipparcos_dat = files('skyshell.catalogues').joinpath('hip_main.dat')
constellation_dat = files('skyshell.catalogues').joinpath('constellation_lines_iau.dat')

def get_sky_plot(
    observation_time    : None | Time =None,
    central_longitude   : float = 0.0, 
    central_latitude    : float = 0.0,
    *,
    plot_hipparcos      : bool = True,
    plot_constellations : bool = True,
    plot_ecliptic       : bool = True,
    plot_equator        : bool = True,
    plot_sun            : bool = True,
    plot_moon           : bool = True,
    plot_daylight       : bool = True,
    plot_moonlight      : bool = False,
    plot_horizon        : bool = False,  
    projection          : str  = 'eckert',
    azimuth             : float= 0.0,
    location            : EarthLocation = PARANAL       
    ):
    """Returns a matplotlib figure of the night sky in galactic coordinates

    Args:
        observation_time (None | Time, optional): Time of observation. Defaults to None.
        central_longitude (float, optional): Galactic Longitude at center of plot. Defaults to 0.0.
        central_latitude (float, optional): Galactic Latitude at center of plot. Only used if projection is 'ortographic'. Defaults to 0.0.
        plot_hipparcos (bool, optional): It True, displays the Hipparcos catalogue. Defaults to True.
        plot_constellations (bool, optional): If True, displays the IAU constellation lines. Defaults to True.
        plot_ecliptic (bool, optional): If True, displays the ecliptic. Defaults to True.
        plot_equator (bool, optional): If True, displays the celestial equator. Defaults to True.
        plot_sun (bool, optional): If True, displays the sun's position. Defaults to True.
        plot_moon (bool, optional): If True, displays the moon's position. Defaults to True.
        plot_daylight (bool, optional): If True, displays simplified daylight contours. Defaults to True.
        plot_moonlight (bool, optional): If True, highlights regions less than 20 degrees from the moon. Defaults to False.
        plot_horizon (bool, optional): If True, displays the horizon line at time of observation. Defaults to False.
        projection (str, optional): Projection to use. Allowed values are ['mercator', 'miller', 'mollweide', 'ortographic', 'eckert']. Defaults to 'eckert'.
        azimuth (float, optional): Azimuth angle. Only used if projection is 'ortographic'. Defaults to 0.0.
        location (EarthLocation, optional): Observation location on earth. Defaults to the Paranal observatory.

    Returns:
        fig, ax: Matplotlib figure and axis 
    """

    #Fetch time
    if observation_time is None:
        time = Time.now()
    else:
        time = observation_time
     
    if not isinstance(time,Time):
        raise ValueError('time must be of class Time')
    
    #Fetch requested projection
    match projection:
        
        case 'mercator':    
            proj = crs.Mercator(
                central_longitude=central_longitude, 
                min_latitude=-80.0, 
                max_latitude=80.0, 
                )

        case 'miller':
            proj = crs.Miller(
                central_longitude=central_longitude, 
                globe=None
                )
        
        case 'mollweide':
            proj = crs.Mollweide(
                central_longitude=central_longitude
                )
             
        case 'ortographic':
            proj = crs.Orthographic(
                central_longitude=central_longitude, 
                central_latitude=central_latitude, 
                azimuth=azimuth)
        case 'eckert':
            proj = crs.EckertIII(
                central_longitude=central_longitude, 
                )        
        case _:
            raise ValueError('Unkown projection')
    
    # Create figure
    fig = plt.figure(figsize=(8, 8), dpi=500)
    ax = plt.axes(projection=proj)

    #Fetch catalogue data
    hipparcos_dict, ra, dec, m = load_hipparcos_positions(hipparcos_dat)
    constellation_lines = load_constellation_polylines(constellation_dat)

    if plot_hipparcos:
             
        #Convert magnitude to flux  
        F0 = 10 
        flux =  F0* 10.0**(-0.4 * m)

        coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")

        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        s =  [min(f,F0/5) for f in flux]
  
        ax.scatter(l, b, s=s, 
                   edgecolors=None, 
                   lw=0, 
                   c='white', 
                   transform=crs.PlateCarree())

    if plot_constellations:
        
        for constellation in constellation_lines:
            for vertex_line in constellation_lines[constellation]:
            
                ra  = [hipparcos_dict[idx][0] for idx in vertex_line]
                dec = [hipparcos_dict[idx][1] for idx in vertex_line]
                
                coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame="icrs")
         
                l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
                b =  coords.galactic.b.deg

                segments = split_on_wrap(l, b, threshold=90)

                for lseg, bseg in segments:
                    ax.plot(lseg, bseg, lw=0.1, c="white", transform=crs.PlateCarree())
            
    if plot_ecliptic:
        
        lambda_ecl = np.linspace(0, 360, 1000) * u.deg
        beta_ecl = np.zeros_like(lambda_ecl)

        ecl_coords = SkyCoord(
            lon=lambda_ecl,
            lat=beta_ecl,
            frame=GeocentricTrueEcliptic
        )

        coords = ecl_coords.transform_to("galactic")
        
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
     
        segments = split_on_wrap(l, b)

        for lseg, bseg in segments:
            ax.plot(lseg, bseg, color="#476EAE", lw=0.8, alpha=0.8, transform=crs.PlateCarree(), zorder=10)

    if plot_equator:

        ra_eq = np.linspace(0, 360, 1000) * u.deg
        dec_eq = np.zeros_like(ra_eq)

        eq_coords = SkyCoord(
            ra=ra_eq,
            dec=dec_eq,
            frame="icrs"
        )

        coords = eq_coords.transform_to("galactic")
 
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
       
        segments = split_on_wrap(l, b)

        for lseg, bseg in segments:
            ax.plot(
                lseg, bseg,
                color="#D97706", 
                lw=0.8,
                alpha=0.8,
                transform=crs.PlateCarree(),
                zorder=10
            )

    if plot_sun:
        
        sun = get_sun_position(time)
          
        coords = sun.transform_to("galactic")
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        
        ax.scatter(l, b, s=60, 
            edgecolors='white', 
            lw=0, 
            c='orange', 
            zorder=10,
            transform=crs.PlateCarree())

    if plot_moon:
            
        moon = get_moon_position(time)

        coords = moon.transform_to("galactic")
        
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        
        ax.scatter(l, b, s=60, 
            edgecolors='white', 
            lw=0, 
            c='white', 
            zorder=10,
            transform=crs.PlateCarree())
 
    if plot_daylight:

        #Sky grid
        n_lon, n_lat = 720, 360
        l = np.linspace(-180, 180, n_lon)
        b = np.linspace(-90, 90, n_lat) 
        L, B = np.meshgrid(l, b)  

        sky = SkyCoord(l=L*u.deg, b=B*u.deg, frame='galactic')

        #Sun separation
        sun = get_sun_position(time)
          
        coords = sun.transform_to("galactic")
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        
        sun = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic')

        sep_sun = sky.separation(sun).deg

        #Plot distinct separations
        ranges = [ (0,90) , (90,96), (96,102), (102,108),(108,112)]
        colors = ["#dceaffff", "#88a6d4ff", "#4575bcff", "#213b66ff", "#192029FF" ]
        angle_transition = 1

        for angle, color in zip(ranges, colors):

            lower_mask = 1 / (1 + np.exp(-(sep_sun - angle[0])/angle_transition))
            upper_mask = 1 / (1 + np.exp(-(sep_sun - angle[1])/angle_transition))

            mask_sun = lower_mask*(1-upper_mask)

            ax.contourf(
                L, B, mask_sun,
                levels=[0.5, 1],
                colors=color,
                alpha=0.8,
                transform=crs.PlateCarree()
            )

    if plot_moonlight:

        #Sky grid
        n_lon, n_lat = 720, 360
        l = np.linspace(-180, 180, n_lon)
        b = np.linspace(-90, 90, n_lat) 
        L, B = np.meshgrid(l, b)  

        sky = SkyCoord(l=L*u.deg, b=B*u.deg, frame='galactic')

        #Moon separation
        moon = get_moon_position(time)
        
        coords = moon.transform_to("galactic")
        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        
        moon = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic')

        sep_moon = sky.separation(moon).deg

        #Plot distinct separations
        ranges = [ (0,10) , (10,15), (15,20)]
        colors = [ "#4575bcff", "#213b66ff", "#192029FF" ]
        angle_transition = 1

        for angle, color in zip(ranges, colors):

            lower_mask = 1 / (1 + np.exp(-(sep_moon - angle[0])/angle_transition))
            upper_mask = 1 / (1 + np.exp(-(sep_moon - angle[1])/angle_transition))

            mask_moon = lower_mask*(1-upper_mask)

            ax.contourf(
                L, B, mask_moon,
                levels=[0.5, 1],
                colors=color,
                alpha=0.8,
                transform=crs.PlateCarree()
            )

    if plot_horizon:
        
        #Sky grid
        n_lon, n_lat = 720, 360
        l = np.linspace(-180, 180, n_lon)
        b = np.linspace(-90, 90, n_lat) 
        L, B = np.meshgrid(l, b)  

        sky = SkyCoord(l=L*u.deg, b=B*u.deg, frame='galactic')

        #Zenith point
        altaz_frame = AltAz(obstime=time, location=location)
        zenith = SkyCoord(alt=90*u.deg, az=0*u.deg, frame=altaz_frame)

        coords = zenith.transform_to("galactic")

        l =  -coords.galactic.l.wrap_at(180 * u.deg).deg
        b =  coords.galactic.b.deg
        
        zenith = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic')

        sep_zenith = sky.separation(zenith).deg

        #Plot distinct separations
        ranges = [ (90,180) ]
        colors = [ "#b75c36ff"]
        angle_transition = 1

        for angle, color in zip(ranges, colors):

            lower_mask = 1 / (1 + np.exp(-(sep_zenith - angle[0])/angle_transition))
            upper_mask = 1 / (1 + np.exp(-(sep_zenith - angle[1])/angle_transition))

            mask_horizon = lower_mask*(1-upper_mask)

            ax.contourf(
                L, B, mask_horizon,
                levels=[0.5, 1],
                colors=color,
                alpha=0.8,
                transform=crs.PlateCarree()
            )

        # Cardinal directions on the horizon (AltAz)
        az_cardinal = [0, 90, 180, 270] * u.deg
        alt_cardinal = [0, 0, 0, 0] * u.deg
        labels = ["N", "E", "S", "W"]

        horizon_cardinal_altaz = SkyCoord(
            az=az_cardinal,
            alt=alt_cardinal,
            frame=altaz_frame
        )

        horizon_cardinal_gal = horizon_cardinal_altaz.transform_to("galactic")

        l_card = -horizon_cardinal_gal.l.wrap_at(180*u.deg).deg
        b_card = horizon_cardinal_gal.b.deg

        ax.scatter(
            l_card,
            b_card,
            s=35,
            c="#b75c36ff",
            zorder=30,
            alpha=0.8,
            transform=crs.PlateCarree()
        )

        for lc, bc, lab in zip(l_card, b_card, labels):
            ax.text(
                 lc, bc,
                 lab,
                 fontsize=4.5,
                 ha="center",
                 va="center",
                 color="black",
                 zorder=31,
                 transform=crs.PlateCarree()
             )

    return fig, ax

