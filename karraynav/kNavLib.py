#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArrayNav.py
Beschreibung: Some standard equations based on WGS-84 model of the earth.
Autor: Luciano Auguto Kruk
Erstellt am: 08.10.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
from math import sqrt, sin
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

    def gravity(self, lat_rad, h_m):
        """
        Using the geografic reference frame, this model allows:
            gravity_n = [0,0,output]
        (see Farrell, (6-141) and (6-142))

        : parameter : lat_rad [rad]  latitude
        : parameter : h_m     [m]    altitude above sea level
        : return    : gravity [m/s2] model of gravity for WGS-84
        """
        s2  = (sin(lat_rad))**2
        s22 = (sin(2.0*lat_rad))**2
        gamma_lat   = 9.780327 * ( 1. + (0.0053024*s2) - (0.0000058*s22) )
        gamma_lat_h = gamma_lat - ((3.0877e-6 - (0.0044e-6*s2) )*h_m) + (0.072e-12*(h_m**2))

        return gamma_lat_h

    def dLat_dt(self, vN, lat_rad, h_m):
        """
        Calculates the derivative of the geografic latitude.

        : parameter : vN      : [m/s] velocity-north
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        : return    : d(latitude)/dt : [rad/s]
        """

        return vN / (self.Rlambda(lat_rad) + h_m)

    def dLong_dt(self, vE, lat_rad, h_m):
        """
        Calculates the derivative of the longitude.

        : parameter : vE      : [m/s] velocity-east
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        : return    : d(longitude)/dt : [rad/s]
        """

        return vE / (cos(lat_rad) * (self.Rphi(lat_rad) + h_m))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
