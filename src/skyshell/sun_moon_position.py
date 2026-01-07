from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import get_body_barycentric

def get_sun_position(time : Time):
    """ Returns the sun's sky coordinates as seen from Earth
    on the specified date 
    """

    earth_pos = get_body_barycentric("earth", time)
    sun_pos   = get_body_barycentric("sun", time)

    sun_from_earth = sun_pos - earth_pos

    sun = SkyCoord(
        sun_from_earth,
        frame="icrs",
        representation_type="cartesian"
    )

    return sun

def get_moon_position(time):
    """ Returns the moon's sky coordinates as seen from Earth
    on the specified date 
    """

    earth_pos = get_body_barycentric("earth", time)
    moon_pos   = get_body_barycentric("moon", time)

    sun_from_earth = moon_pos - earth_pos

    moon = SkyCoord(
        sun_from_earth,
        frame="icrs",
        representation_type="cartesian"
    )

    return moon
