import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from skyshell.skyplot import get_sky_plot

time = Time("2026-05-28 06:00:00", scale="utc")
location = EarthLocation(lat=-24.6272*u.deg, lon=-70.4042*u.deg, height=2635*u.m)

#For the current time
fig, ax = get_sky_plot(
    observation_time    = time,
    central_longitude   = 0.0, 
    central_latitude    = 0.0,
    plot_hipparcos      = True,
    plot_constellations = True,
    plot_ecliptic       = True,
    plot_equator        = True,
    plot_sun            = True,
    plot_moon           = True,
    plot_daylight       = True,
    plot_moonlight      = False,
    plot_horizon        = False,  
    projection          = 'eckert',
    azimuth             = 0.0,
    location            = location  
    )
    
fig.savefig(f'output.png', bbox_inches='tight')