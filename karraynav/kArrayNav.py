#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArrayNav.py
Beschreibung: Mathematical manipulation of vectors and matrices as needed by navigation modelling, and more.
Autor: Luciano Auguto Kruk
Erstellt am: 06.09.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
from math import sin, cos, tan, sqrt, pi, atan, atan2, asin
from mpmath import sec
from .kArray import kArray
from .kNavLib import kNavLib
#import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayLib:
    def to_skew(self):
        vtype = self._type(self.array)
        if vtype in [ self.TYPE_HORIZONTAL, self.TYPE_VERTICAL ]:
            val  = self.array.squeeze()
            assert len(val) == 3
            return self.__class__( [
                [0, -val[2], val[1]],
                [val[2], 0, -val[0]],
                [-val[1], val[0], 0] 
            ])
        else:
            raise(NameError("this is not a valid vector 3x1"))

    def X(self, y):
        """
        Calculates the cross-product between self and y.
        : input  : y = [3x1] vector
        : return : self X y
        """
        assert self.shape == (3,1)
        assert y.shape == (3,1)
        a = self.array.squeeze()
        b = y.array.squeeze()

        return self.__class__( [
            (a[1]*b[2]) - (a[2]*b[1]),
            (a[2]*b[0]) - (a[0]*b[2]),
            (a[0]*b[1]) - (a[1]*b[0])
        ], hvector=False )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kNavTransformations(kNavLib):

    def to_deg(self):
        return self.__class__(self.array * 180. / pi)
        #return val * 180. / pi

    def to_rad(self):
        return self.__class__(self.array * pi / 180.)
        #return val * pi / 180.

    def euler2Q(self):
        """
        Navigation -- from euler to Q.
        : input: kArray = [phi, theta, psi]
        : parameter : phi   [rad]
        : parameter : theta [rad]
        : parameter : psi   [rad]
        : return: kArray = Q4
        """
        array = self.array.squeeze()
        assert len(array) == 3
        phi   = array[0]
        theta = array[1]
        psi   = array[2]

        half_phi   = 0.5*phi
        half_theta = 0.5*theta
        half_psi   = 0.5*psi;

        return self.__class__( [
            (cos(half_phi)*cos(half_theta)*cos(half_psi)) + (sin(half_phi)*sin(half_theta)*sin(half_psi)),
            (sin(half_phi)*cos(half_theta)*cos(half_psi)) - (cos(half_phi)*sin(half_theta)*sin(half_psi)),
            (cos(half_phi)*sin(half_theta)*cos(half_psi)) + (sin(half_phi)*cos(half_theta)*sin(half_psi)),
            (cos(half_phi)*cos(half_theta)*sin(half_psi)) - (sin(half_phi)*sin(half_theta)*cos(half_psi))
        ], hvector=False );

    def Q2euler(self):
        """
        Navigation -- from Q to euler.
        : input: kArray = Q4
        : output   : phi   [rad]
        : output   : theta [rad]
        : output   : psi   [rad]
        : return: kArray( [phi, theta, psi] )
        """

        q     = self.array.squeeze()
        assert len(q) == 4

        phi   = atan2(2.0*((q[2]*q[3])+(q[0]*q[1])), (q[0]**2.0)-(q[1]**2.0)-(q[2]**2.0)+(q[3]**2.0));
        psi   = atan2(2.0*((q[1]*q[2])+(q[0]*q[3])), (q[0]**2.0)+(q[1]**2.0)-(q[2]**2.0)-(q[3]**2.0));
        try:
            theta = asin(2.0*((q[0]*q[2])-(q[1]*q[3])));
        except ValueError:
            print("ERR: norm(Q) = {:f}".format(np.sqrt(np.sum(q**2))))
            theta = 0;

        return self.__class__( (phi, theta, psi) )

    def Q2C(self):
        """
        Navigation -- from Q to C.

        If Q represents the transformation from 'a' to 'b', the matrix
        'C' represents 'Ca2b'.

        : input    : kArray = q
        : output   : C
        """

        q = self.array.squeeze()
        assert len(q) == 4

        C = kArray( np.empty((3,3)) )
        C[0,0] = (q[0]**2.0) + (q[1]**2.0) - (q[2]**2.0) - (q[3]**2.0);
        C[0,1] = 2.0 * ((q[1]*q[2]) + (q[0]*q[3]));
        C[0,2] = 2.0 * ((q[1]*q[3]) - (q[0]*q[2]));

        C[1,0] = 2.0 * ((q[1]*q[2]) - (q[0]*q[3]));
        C[1,1] = (q[0]**2.0) - (q[1]**2.0) + (q[2]**2.0) - (q[3]**2.0);
        C[1,2] = 2.0 * ((q[2]*q[3]) + (q[0]*q[1]));

        C[2,0] = 2.0 * ((q[1]*q[3]) + (q[0]*q[2]));
        C[2,1] = 2.0 * ((q[2]*q[3]) - (q[0]*q[1]));
        C[2,2] = (q[0]**2.0) - (q[1]**2.0) - (q[2]**2.0) + (q[3]**2.0);

        return self.__class__( C )

    def C2euler(self):
        """
        Navigation -- from C to (phi,theta,psi)[rad]
        """

        C = self.array
        assert C.shape == (3,3)
        assert(C[2,2] != 0)
        assert(C[0,0] != 0)
        assert(C[0,2]>=-1 and C[0,2]<=1)

        phi   = np.arctan2(C[1,2], C[2,2])
        theta = np.arcsin(-C[0,2])
        psi   = np.arctan2(C[0,1], C[0,0])

        return self.__class__( (phi, theta, psi) )

    def euler2C(self):
        """
        : Convert euler [rad] to C matrix.
        : order = [phi, theta, psi] [rad]
        """

        a = self.array.squeeze()
        assert len(a) == 3

        if False:
            phi    = a[0]
            theta  = a[1]
            psi    = a[2]
            cphi   = cos(phi)
            ctheta = cos(theta)
            cpsi   = cos(psi)
            sphi   = sin(phi)
            stheta = sin(theta)
            spsi   = sin(psi)

            C = self.__class__( np.empty((3,3)) )
            C[0][0] = cpsi*ctheta
            C[0][1] = spsi*ctheta
            C[0][2] = -stheta
            C[1][0] = (-spsi*cphi) + (cpsi*stheta*sphi)
            C[1][1] = (cpsi*cphi) + (spsi*stheta*sphi)
            C[1][2] = ctheta * sphi
            C[2][0] = (sphi*sphi) + (cpsi*stheta*cphi)
            C[2][1] = (-cpsi*sphi) + (sphi*stheta*cphi)
            C[2][2] = ctheta * cphi

            #return C

        return self.__class__( a ).euler2Q().Q2C()

    def C2Q(self):
        """
        Navigation -- from C to Q
        """
        return self.__class__( self.C2euler().euler2Q() )

    def Qinv(self):
        """
        Calculates the inverse of a quaternion. This is similar to inverting a transformation matrix.
        """
        shape = self.array.shape
        tmp   = self.array.squeeze().tolist()
        ret   = [tmp[0]] + [-i for i in tmp[1:]]
        return self.__class__( np.asarray(ret).reshape(shape) / sum( [i**2 for i in ret] ))

    def q1_x_q2(self, q2):
        """
        Navigation -- multiplies two quaternions

        Let q1 represent C_a2b, and q2 represent C_b2c.
        The product C_a2c = C_b2c.C_a2b might be represented
        by q_a2c = q_b2c.q_a2b

        output: np.array quaternion q3=q1.q2
        """
        q1 = self.array.squeeze()
        assert len(q1) == 4

        if isinstance(q2, kArrayNav) or isinstance(kArray):
            q2 = q2.array.squeeze()
        assert len(q2) == 4

        q3 = np.array([
            (q1[0]*q2[0])-(q2[1]*q1[1])-(q2[2]*q1[2])-(q2[3]*q1[3]),
            (q2[0]*q1[1])+(q2[1]*q1[0])+(q2[2]*q1[3])-(q2[3]*q1[2]),
            (q2[0]*q1[2])+(q2[2]*q1[0])-(q2[1]*q1[3])+(q2[3]*q1[1]),
            (q2[0]*q1[3])+(q2[3]*q1[0])+(q2[1]*q1[2])-(q2[2]*q1[1])
        ])

        return self.__class__( q3 )

    def Re2n(self, lat, lon):
        """
        Navigation -- calculates Re2n(lat,lon)
        The result (Re2n) does not change the content of the current object.

        : input    : lat   [rad]
        : input    : lon   [rad]
        : output   : Re2n
        """

        Re2n = np.empty((3,3));
        Re2n[0,0] = -sin(lat)*cos(lon);
        Re2n[0,1] = -sin(lat)*sin(lon);
        Re2n[0,2] = cos(lat);
        Re2n[1,0] = -sin(lon);
        Re2n[1,1] = cos(lon);
        Re2n[1,2] = 0;
        Re2n[2,0] = -cos(lat)*cos(lon);
        Re2n[2,1] = -cos(lat)*sin(lon);
        Re2n[2,2] = -sin(lat);

        return self.__class__( Re2n )

    def ecef_llh2xyz(self):
        """
        : Convert from ECEF-geodetic to XYZe.
        : input   : llh  = (lat_rad, lon_rad, h_m)
        : output  : xzy_e [m]
        """

        geo = self.array.squeeze()
        assert len(geo) == 3
        lat = geo[0]
        lon = geo[1]
        h   = geo[2]

        s  = sin(lat);
        RN = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));

        return self.__class__( [
            (RN + h) * cos(lat) * cos(lon),
            (RN + h) * cos(lat) * sin(lon),
            ((RN * (1.0 - self.earth_e2)) + h) * sin(lat)
        ] )

    def ecef_xyz2llh(self):
        """
        : Convert from XYZe to ECEF-geodetic.
        : input  : xyz_e [m]
        : output : llh  = (lat_rad, lon_rad, h_m)
        """

        rect = self.array.squeeze()
        assert len(rect) == 3
        x = rect[0]
        y = rect[1]
        z = rect[2]

        p = sqrt((x * x) + (y * y));
        h     = 0;
        RN          = self.earth_a;
        for i in range(100): # timeout
            lasth   = h;

            # algorithm (Farrell/Barth p.28)
            s   = z / (((1.0 - self.earth_e2) * RN) + h);
            lat = atan((z + (self.earth_e2 * RN * s)) / p);
            RN  = self.earth_a / sqrt(1.0 - (self.earth_e2 * s * s));
            h   = (p / cos(lat)) - RN;

            # centimeter accuracy:
            if abs(lasth-h) < 0.01:
               break;

        lon = atan2(y, x);

        return self.__class__( [lat, lon, h] )

    def gravity_n(self, lat_rad, h_m):
        """
        Calculates the local gravity vector in the geografic frame (n).

        : parameter : lat_rad [rad]  latitude
        : parameter : h_m     [m]    altitude above sea level
        : return    : vector with local gravity [3x1]
        """

        return self.__class__( [0,0,self.gravity(lat_rad, h_m)], hvector=False )

    def dqdt(self, w):
        """
        The derivative of the quaternions is $\dot{q} = 1/2 .B(w).q$
        This funtion returns $\dot{q}$.
        : input  : q4
        : output : d(q4)dt
        """

        K      = 1e1
        cq     = kArray( [i for i in self], hvector=False )
        cq2    = kArray( [i*i for i in self], hvector=False )
        W      = kArray( w, hvector=True )
        epslon = 1.0 - sum(cq2.to_list())

        B = kArray( [
            [   0, -W[0][0], -W[0][1], -W[0][2]],
            [W[0][0],     0,  W[0][2], -W[0][1]],
            [W[0][1], -W[0][2],     0,  W[0][0]],
            [W[0][2],  W[0][1], -W[0][0],     0]
        ] )

        dq = (0.5 * B * cq) + (K*epslon*cq)
        return self.__class__( dq )

    def dEulerDt(self, w):
        """
        Calculates the derivative vector of the euler angles.
        See Titerton (3-52).
        : input  : [phi, theta, psi] [rad]
        : return : d[phi, theta, psi]/dt [rad/s]
        """

        a = self.array.squeeze()
        assert len(a) == 3
        phi    = a[0]
        theta  = a[1]
        psi    = a[2]

        if isinstance(w, list):
            b = w
        elif isinstance(w, self.__class__):
            b = w.array.squeeze().tolist()
        else:
            raise(NameError("this type is still not supported"))

        aux = (b[1]*sin(phi)) + (b[2]*cos(phi))
        dphi = (aux * tan(theta)) + b[0]
        dtta = (b[1]*cos(phi)) - (b[2]*sin(phi))
        dpsi = aux * sec(theta)

        return self.__class__( [dphi, dtta, dpsi] )

    def dLLH_dt(self, vN, vE, vD, lat_rad, h_m):
        """
        Calculates the vector with the derivatives of the
        latitude, longitude and altitude.

        : parameter : vN      : [m/s] velocity-north
        : parameter : vE      : [m/s] velocity-east
        : parameter : vD      : [m/s] velocity-down
        : parameter : lat_rad : [rad] latitude
        : parameter : h_m     : [m]   altitude above sea level
        """

        return self.__class__( [
            self.dLat_dt(vN, lat_rad, h_m),
            self.dLong_dt(vE, lat_rad, h_m),
            -vD ], hvector=False )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayNav (kArray, kArrayLib, kNavTransformations):
    def __init__(self, *args, **kargs):
        if len(args) > 0:
            super().__init__(*args, **kargs)

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayNavTests:

    def do_tests(self):
        import scipy.integrate  as Int

        print("==== to_skew() ====")
        a = kArrayNav( [1,2,3], hvector=False )
        b = [0,-3,2,3,0,-1,-2,1,0]
        assert all([i==j for i,j in zip(a.to_skew(),b)])
        print(a.to_skew())

        print("==== cross product ====")
        a = kArrayNav( [5,-2,4], hvector=True )
        b = kArrayNav( [-1,2,3], hvector=True )

        try:
            ok = False
            c  = a.X(b)
        except:
            ok = True
        if not ok:
            raise(NameError("should not get here!"))

        c1 = a.T.X(b.T)
        c2 = np.cross(a().squeeze(), b().squeeze())
        c3 = a.to_skew() * b.T

        for i,j,k in zip(c1,c2,c3):
            assert abs(i-j) < 1e-10
            assert abs(i-k) < 1e-10

        print("==== to_rad, to_deg ====")
        euler_deg   = kArrayNav( [-30, 10, 170] )
        euler_rad   = euler_deg.to_rad()
        euler_deg_t = euler_rad.to_deg()

        for i,j in zip(euler_deg, euler_deg_t):
            assert abs(i-j) < 1e-10

        print("==== euler - Q - euler ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            euler   = kArrayNav( euler_rad )
            q4      = euler.euler2Q()
            euler_t = q4.Q2euler()
            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

        print("==== Q2C / C2Q ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            euler   = kArrayNav( euler_rad )
            q4    = euler.euler2Q()
            C     = q4.Q2C()
            euler_t = C.C2euler()
            q4_t  = C.C2Q()

            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

            for j,k in zip(q4, q4_t):
                assert abs(j-k) < 1e-10

        print("==== q1_x_q2 ====")
        for i in range(20):
            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            Ca2b = kArrayNav( euler_rad ).euler2C()
            qa2b = kArrayNav( euler_rad).euler2Q()

            euler_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            Cb2c = kArrayNav( euler_rad ).euler2C()
            qb2c = kArrayNav( euler_rad ).euler2Q()

            Ca2c = Cb2c * Ca2b
            qa2c = qb2c.q1_x_q2(qa2b)

            euler   = Ca2c.C2euler()
            euler_t = qa2c.Q2euler()

            for j,k in zip(euler, euler_t):
                assert abs(j-k) < 1e-10

        print("==== Re2n ====")
        Re2n = kArrayNav().Re2n(0,0)
        assert Re2n == kArrayNav( [[0,0,1], [0,1,0], [-1,0,0]] )

        print("==== llh / xyz_e ====")
        for i in range(20):
            angles_rad = (((np.random.rand(1,3) * 2.0) - 1.0) * 90.) * pi/180
            llh   = kArrayNav( angles_rad )
            xyz   = llh.ecef_llh2xyz()
            llh_t = xyz.ecef_xyz2llh()
        
            for i,j in zip(llh, llh_t):
                #print("{:f} == {:f} ?".format(i,j))
                assert abs(i-j) < 1e-3

        #----------------------#
        # some dynamic tests:
        #----------------------#
        print("==== dqdt ====")
        from   scipy.integrate import odeint;
        from   numpy           import dot;
        print()

        #  I: inertial frame
        #  b: body frame
        qI2b = kArrayNav( [0,0,0] ).euler2Q()

        # angular rotation between I and b:
        # \omega_{Ib}^I
        w_ib_i = kArray( [2.*pi/180., 0, 0] )

        def eqdiff(q,t,w_ib_i):
            """
            we will use scipy to call this function, and therefore q shall be a list().
            """
            qi2b = kArrayNav(q)
            RI2b = qi2b.Q2C()
            dqdt = qi2b.dqdt( RI2b * kArray( w_ib_i, hvector=False ))
            return dqdt.to_list()

        # a vector described at I:
        F_i = kArray( [0,0,1], hvector=False )
        print("F_i = ")
        print(F_i)

        for t in [1,5,20,90]:
            print()
            # after t seconds, the quaternions should be:
            #print("qI2b.to_list() =")
            #print(qI2b.to_list())
            y = odeint(eqdiff, qI2b.to_list(), [0,t], (w_ib_i,))[1,:]
            # with these euler angles:
            euler = ((180./pi) * kArrayNav(y).Q2euler()).to_list()
            print('euler = {:s}'.format(str(euler)))

            # and described at b:
            F_b = (kArrayNav(y).Q2C() * F_i).to_list()

            print("F_b(phi = {:1.03f}) = [{:1.03f} {:1.03f} {:1.03f}]".format(
                euler[0], F_b[0], F_b[1], F_b[2]))

        # gravity:
        import matplotlib.pylab as plt

        plt.figure(1).clf()
        fig, ax = plt.subplots(1,1,num=1)
        tmp = kArrayNav([0])
        leg = list()
        for lat in kArrayNav( [0, 45, 80] ).to_rad():
            leg.append("{:1.1f}".format(lat*180/pi))
            g = list()
            for h in np.linspace(0, 3000, 20):
                g.append( tmp.gravity(lat, h) )
            ax.plot(np.linspace(0,3000,20), g)

        ax.legend(leg)

        #----------------------#
        # coherence tests
        # (sense among frames)
        #----------------------#
        list_euler = [
                kArrayNav( [0,0,45] ).to_rad(),
                kArrayNav( [0,45,0] ).to_rad(),
                kArrayNav( [45,0,0] ).to_rad(),
        ]

        # with the euler angles above, the transformation of vector [1,0,0]_b will lead to these results:
        sq2 = sqrt(2.)/2.
        list_results_at_n = [
                kArrayNav( [sq2, sq2, 0], hvector=False ),
                kArrayNav( [sq2, 0, -sq2], hvector=False ),
                kArrayNav( [1,0,0], hvector=False ),
        ]

        # Tests #1
        for euler,vn in zip(list_euler, list_results_at_n):
            # from a fixed reference frame 'n' to 'b'
            Cn2b = euler.euler2C()

            # transforming a vector from 'b' to 'n':
            rb = kArrayNav( [1,0,0], hvector=False )
            rn = Cn2b.T * rb
            assert rn == vn

        # Tests #2
        # An angular velocity [0,0,1] shall increase \psi
        # The angle \psi refers to how 'n' sees 'b'.
        Cn2b   = kArrayNav( [0,0,0] ).to_rad().euler2C()
        Cb2n   = Cn2b.T # to use w_nb_b

        w_nb_b = kArrayNav( [0,0,1], hvector=False )
        Cb2n_p = Cb2n * w_nb_b.to_skew()   # derivative of a tranformation matrix
        Cb2n_step = Cb2n + (0.01 * Cb2n_p) # small integration step
        euler  = Cb2n_step.T.C2euler()     # transposing to obtain the euler again
        assert euler[0][2] > 0
        
        # Tests #3
        # An angular velocity [1,1,1] shall increase all euler angles.
        # The euler angles refer to how 'n' sees 'b'.
        Cn2b   = kArrayNav( [0,0,0] ).to_rad().euler2C()
        Cb2n   = Cn2b.T # to use w_nb_b

        w_nb_b = kArrayNav( [1,1,1], hvector=False )
        Cb2n_p = Cb2n * w_nb_b.to_skew()   # derivative of a tranformation matrix
        Cb2n_step = Cb2n + (0.01 * Cb2n_p) # small integration step
        euler  = Cb2n_step.T.C2euler()     # transposing to obtain the euler again
        for i in euler[0]:
            assert i > 0

        # Tests #4
        # Similar to above, but with quaternions in parallel.
        # Note: w_nb is described at b.
        w_nb_b = kArrayNav( [5,5,5], hvector=False ).to_rad() # [rad/s]

        # initial euler angles, from n to b:
        euler_n2b = kArrayNav( [10,10,10] ).to_rad()

        # [k]
        Rb2n = euler_n2b.euler2C().T
        qn2b = euler_n2b.euler2Q()

        # d../dt:
        Rb2n_p = Rb2n * w_nb_b.to_skew() # to use w_nb_b, here the derivative needs Rb2n
        qn2b_p = qn2b.dqdt(w_nb_b)       # to use w_nb_b, here the derivative needs qn2b

        # one small integration step [k+1] (first order):
        dt = 0.01 # [s]
        Rb2n += dt * Rb2n_p
        qn2b += dt * qn2b_p

        # euler:
        euler_k1_R_b2n = Rb2n.T.C2euler().to_deg() # transposing to get Rb2n
        euler_k1_q_n2b = qn2b.Q2euler().to_deg()

        # tests:
        for i in range(3):
            assert euler_k1_R_b2n[0][i] > euler_n2b[0][i]
            assert euler_k1_q_n2b[0][i] > euler_n2b[0][i]
            assert abs( euler_k1_R_b2n[0][i] - euler_k1_q_n2b[0][i] ) < 1e-4

        # Tests #5
        # Here we want to check whether inverting a quaternion has the same effect as
        # inverting a transformation matrix.
        for i in range(20):
            phi   = 20*np.random.randn()
            theta = 20*np.random.randn()
            psi   = 20*np.random.randn()

            euler = kArrayNav( [phi, theta, psi] ).to_rad()
            R = euler.euler2C()

            q = euler.euler2Q()
            q_inv = q.Qinv()

            euler_from_q_inv = q_inv.Q2euler().to_deg()
            euler_from_R_inv = R.T.C2euler().to_deg()
            assert euler_from_q_inv == euler_from_R_inv

            if False:
                print("euler angles at [k]")
                print(euler.to_deg())
                print("euler angles at [k+1] (from q)")
                print(euler_from_q_inv)
                print("euler angles at [k+1] (from R)")
                print(euler_from_R_inv)
                print()

        # Tests #6
        # To compare the effect of w_nb_b on the euler angles using two methods:
        #  1) dEulerDt()
        #  2) Rb2n_p = Rb2n * w_nb_b
        phi   = 20*np.random.randn()
        theta = 20*np.random.randn()
        psi   = 20*np.random.randn()

        def test6_eqdiff_dEulerDt(y, t, *args):
            w = args[0](t)
            dEuler = kArrayNav(y).dEulerDt(w).to_list()
            return dEuler

        def test6_eqdiff_dRdt(y, t, *args):
            w = args[0](t)
            R = kArrayNav( np.asarray(y).reshape(3,3) )
            Rp = R * kArrayNav(w).to_skew()
            return Rp.to_list()

        def tests_angular_speed(t):
            if t < 1:
                w   = kArrayNav( [10,10,10] ).to_rad().to_list()
            elif t < 2:
                w   = kArrayNav( [-10,10,-10] ).to_rad().to_list()
            else:
                w   = kArrayNav( [0,-10,0] ).to_rad().to_list()
            return w

        T   = np.linspace(0, 5, 5*100)

        ret = Int.odeint( test6_eqdiff_dEulerDt, kArrayNav([phi, theta, psi]).to_rad().to_list(), T, args=(tests_angular_speed,) )
        ret = [kArrayNav(i).to_deg().to_list() for i in ret]

        Ret = Int.odeint( test6_eqdiff_dRdt, kArrayNav([phi, theta, psi]).to_rad().euler2C().T.to_list(), T, args=(tests_angular_speed,) )
        Ret = [(kArrayNav( np.asarray(i).reshape(3,3) ).T.C2euler().to_deg() + 0.1).to_list() for i in Ret]
        # this 0.1 above is just to separate the curves on the next graph

        plt.figure(2).clf()
        fig,ax = plt.subplots(1,1,num=2)
        ax.plot(T, ret)
        ax.plot(T, Ret)
        ax.grid(True)

        for i,j in zip(ret, Ret):
            for m,n in zip(i,j):
                assert abs(m-n) < 0.11

        ###########################
        plt.show(block=False)
        ###########################
        print("ok")

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
if __name__ == "__main__":
    tests = kArrayNavTests()
    tests.do_tests()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
