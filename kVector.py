#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kVector.py
Beschreibung: Mathematical manipulation of vectors as needed by navigation modelling.
Autor: Luciano Auguto Kruk
Erstellt am: 27.08.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
import math
import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kVector:
    def __init__(self, val, transposed=False):
        # by default, neveer trasposed:
        self.transposed = transposed

        # loading:
        if isinstance(val, kVector):
            self.vector = copy.deepcopy(val.vector)
        elif isinstance(val, list) or isinstance(val, tuple):
            self.vector = list(copy.deepcopy(val))
        else:
            self.vector = None
            raise("format of 'val' is unknown")

    def _do_format(self, fmt):
        txt = "[ {{:{:s}}}; {{:{:s}}}; {{:{:s}}} ]".format(fmt, fmt, fmt).format( self.vector[0], self.vector[1], self.vector[2] )
        if self.transposed:
            return txt + ".T"
        else:
            return txt

    def __repr__(self):
        return "<class kVector {:s}>".format(self._do_format("f"))

    def __format__(self, fmt):
        return self._do_format(fmt)

    @property
    def T(self):
        return kVector(self.vector, transposed=True)

    #( --- sum --- )#
    def __add__(self, y):
        if isinstance(y, kVector):
            ret = kVector([ self.vector[0] + y.vector[0], self.vector[1] + y.vector[1], self.vector[2] + y.vector[2] ])
        elif isinstance(y, int) or isinstance(y, float):
            ret = kVector([ self.vector[0] + y, self.vector[1] + y, self.vector[2] + y ])
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))
        return ret

    def __radd__(self, y):
        return self.__add__(y)

    def __iadd__(self, y): # +=
        q       = self.__add__(y)
        self.vector = q.vector
        return(self)

    #( --- negative signal --- )#
    def __neg__(self):
        return kVector([-i for i in self.vector])

    #( --- difference --- )#
    def __sub__(self, y):
        if isinstance(y, kVector):
            ret = kVector([ self.vector[0] - y.vector[0], self.vector[1] - y.vector[1], self.vector[2] - y.vector[2] ])
        elif isinstance(y, int) or isinstance(y, float):
            ret = kVector([ self.vector[0] - y, self.vector[1] - y, self.vector[2] - y ])
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))
        return ret

    def __rsub__(self, y):
        return self.__add__(-y)

    def __isub__(self, y): # -=
        q = self.__sub__(y)
        self.vector = q.vector
        return(self)

    #( --- multiplication --- )#
    def __mul__(self, y):
        if isinstance(y, kVector):
            if self.transposed:
                # [1x3] x [3x1]
                ret = (self.vector[0] * y.vector[0]) + (self.vector[1] * y.vector[1]) + (self.vector[2] * y.vector[2])
            elif y.transposed:
                # [3x1] x [1x3]
                raise(NameError("a matrix class needs to be developed"))
            else:
                raise(NameError("WTF??"))

        elif isinstance(y, int) or isinstance(y, float):
            ret = kVector([ y*self.vector[0], y*self.vector[1], y*self.vector[2] ])

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

    def __rmul__(self, y):
        return self.__mul__(y)

    def __imul__(self, y): # *=
        q = self.__mul__(y)
        self.vector = q.vector
        return(self)

    #( --- division --- )#
    def __truediv__(self, y):
        raise(NameError("WTF??"))

    def __rtruediv__(self, y):
        raise(NameError("WTF??"))

    def __itruediv__(self, y): # /=
        raise(NameError("WTF??"))

    #( --- miscelaneous --- )#
    def __abs__(self):
        return math.sqrt( (self.vector[0]**2) + (self.vector[1]**2) + (self.vector[2]**2) )

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
if __name__ == "__main__":
    print("==== __init__() ====")
    a = kVector((1,2,3))
    print(a)
    print("a = {:f}".format(a))

    b = kVector(a)
    print(b)
    print("b = {:f}".format(b))

    print("==== sum ====")
    a = kVector((1,2,3))
    b = kVector((4,5,6))
    print("a + b    = {:f}".format(a+b))
    print("a + (-1) = {:f}".format(a+(-1)))
    print("(-1) + a = {:f}".format((-1)+a))
    a += 7
    print("a += 7: {:f}".format(a))

    print("==== negative signal ====")
    print("-a = {:f}".format(-a))

    print("==== difference ====")
    a = kVector((1,2,3))
    b = kVector((4,5,6))
    print("a - b   = {:f}".format(a-b))
    print("3.0 - a = {:f}".format(3.0-a))
    a -= 3.0
    print("a -= 3.0: {:f}".format(a))

    print("==== multiplication ====")
    a = kVector([1,2,3])
    b = kVector((4,5,6))
    try:
        print("a * b = {:f}".format(a*b))
    except:
        print("a * b =")
        print("    ^^^ERROR: 'a' should be transposed!")

    print("a.T     = {:f}".format(a.T))
    print("a.T * b = {:f}".format( a.T*b ))
    print("7.0 * a = {:f}".format( 7.0 * a ))
    a *= 10
    print("a *= 10 : {:f}".format(a))

    print("==== norm ====")
    a = kVector([1,2,3])
    print("||a|| = {:f}".format(abs(a)))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
