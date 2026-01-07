import ast 
import numpy as np

def load_constellation_polylines(filename):
    """Returns a dictionary with the constelation's asterisms as:
    
        {
            "Andromeda": [
                [hip1, hip2, hip3, ...],
                [hipA, hipB, ...],
            ],
            ...
        }
    
    The constelations are defined according to the IAU convention and the
    file used corresponds to the one obtained from:

    (https://github.com/dcf21/constellation-stick-figures)
    
    """
    constellations = {}
    current_constellation = None

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
        
            if not line or line.startswith("#"):
                continue

            # Constellation header
            if line.startswith("*"):
                current_constellation = line[1:].strip()
                constellations[current_constellation] = []
                continue

            # Polyline definition
            if line.startswith("[") and current_constellation is not None:
                hip_list = ast.literal_eval(line)
                hip_list = [int(h) for h in hip_list]
                constellations[current_constellation].append(hip_list)

    return constellations

def load_hipparcos_positions(filename):
    """Returns a dictionary with the Hipparcos catalogue stars
    and corresponding (RA,DEC) coordinates. The dictionary has
    the format:

   { HIP : (RA_deg, Dec_deg) }

    """
    hip_pos = {}
    
    ra_list = []
    dec_list = []
    mag_list = []


    with open(filename, "r") as f:
        for line in f:
            fields = [f.strip() for f in line.split("|")]

            hip = int(fields[1])
            ra  = fields[3]
            dec = fields[4]
            mag = fields[5]

            try:
                mag = float(mag)
            except:
                continue
            
            h,m,s = [float(val) for val in ra.split()]
            ra = 15 * (h + m/60 + s/3600)

            # print(dec.split())
            deg,m,s = [float(val) for val in dec.split()]
            dec = deg + m/60 + s/3600

            ra_list.append(ra)
            dec_list.append(dec)
            mag_list.append(mag)
            
            hip_pos[hip] = (ra, dec, mag)

    return hip_pos, np.array(ra_list), np.array(dec_list), np.array(mag_list)

def split_on_wrap(lrad, brad, threshold=180):
    """Given a list of points, returns the unwrapped
    values that allow proper plotting on the sky plot
    """
    segments = []
    start = 0
    for i in range(1, len(lrad)):
        if abs(lrad[i] - lrad[i-1]) > threshold:
            segments.append((lrad[start:i], brad[start:i]))
            start = i

    segments.append((lrad[start:], brad[start:]))
    return segments
