#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kNavLib.py
Beschreibung: Library with some relevant transformations for navigation math.
Autor: Luciano Auguto Kruk
Erstellt am: 09.09.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
from math import atan, atan2, sin, cos, sqrt
#from kArray import kArray
#import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kNavLib:
    # Earth Elliptic Model #
    earth_a  = 6378137.0; # [m]
    earth_b  = 6356752.3142; # [m]
    wie      = 1.0 * 7.2921151467e-5;
    earth_f  = (earth_a-earth_b)/earth_a;
    earth_e  = sqrt(earth_f*(2.0-earth_f));
    earth_e2 = (earth_e**2.0);

    def Rlambda(self, lat_rad):
        """
        : parameter : lat_rad [rad] latitude
        : output    : R_lbd
        : (Farrell/Barth, eq 6-13)
        """
        return (self.earth_a*(1.-self.earth_e2)) / ((1.-(self.earth_e2*(sin(lat_rad)**2)))**1.5);

    def Rphi(self, lat_rad):
        """
        : parameter : lat_rad [rad] latitude
        : output    : R_phi
        : (Farrell/Barth, eq 6-14)
        """
        return self.earth_a / sqrt(1.-(self.earth_e2*(sin(lat_rad)**2.0)));

    def llh2xyze(self, llh):
        """
        : Convert from ECEF-geodetic to XYZe.
        : parameter : llh  = (lat_rad, lon_rad, h_m)
        : return    : XYZe = (x_m, y_m, z_m)
        """

        lat, lon, h = llh
        s  = sin(lat);
        RN = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));

        x = (RN + h) * cos(lat) * cos(lon)
        y = (RN + h) * cos(lat) * sin(lon)
        z = ((RN * (1.0 - self.earth_e2)) + h) * sin(lat)

        return (x,y,z)

    def xzye2llh(self, xyz_e):
        """
        : Convert from XYZe to ECEF-geodetic.
        : parameter : xyz_e [m]
        : return    : llh  = (lat_rad, lon_rad, h_m)
        """

        x,y,z = xyz_e
        h     = 0
        RN    = self.earth_a
        p     = sqrt((x * x) + (y * y))
        for i in range(100): # timeout
            #print "[lat h] = [%1.09f %1.03f]" % (lat, h)
            lasth = h;

            # algorithm (Farrell/Barth p.28)
            s   = z / (((1.0 - self.earth_e2) * RN) + h);
            lat = atan((z + (self.earth_e2 * RN * s)) / p);
            RN  = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));
            h   = (p / cos(lat)) - RN;

            # centimeter accuracy:
            if abs(lasth-h) < 0.01:
               break;

        lon = atan2(y, x);

        return lat, lon, h

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
if __name__ == "__main__":

    print("==== llh / xyz_e ====")
    nav = kNavLib()
    lat, lon, h = nav.xzye2llh( nav.llh2xyze( (10./57., -30./57., 33) ))
    assert abs(lat - 10./57) < 1e-5
    assert abs(lon + 30./57) < 1e-5
    assert abs(h   - 33.   ) < 1e-2
    print("lat = {:f}".format(lat*57.))
    print("lon = {:f}".format(lon*57.))
    print("h   = {:f}".format(h))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
