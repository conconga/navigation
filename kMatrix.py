#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kMatrix.py
Beschreibung: Mathematical manipulation of matrices as needed by navigation modelling.
Autor: Luciano Auguto Kruk
Erstellt am: 27.08.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kMatrix:
    """
    The matrix is stored as a vector with the elements fill rows before columns.
    """
    def __init__(self, val, transposed=False, size=(3,3)):
        # by default, never trasposed:
        self.transposed = transposed

        # by default, [3x3]:
        self.size = size

        # loading:
        if isinstance(val, kMatrix):
            self.size       = copy.deepcopy(val.size)
            self.vector     = copy.deepcopy(val.vector)

        elif isinstance(val, list) or isinstance(val, tuple):
            self.vector    = list(copy.deepcopy(val))

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        # nb of items:
        self.len  = size[0] * size[1]

    def __getitem__(self, index):
        """
        index = (row, column)
        """
        if self.transposed:
            idx = (index[1] * self.size[1]) + index[0]
        else:
            idx = (index[0] * self.size[1]) + index[1]

        return self.vector[idx]

    def _do_format(self, fmt):
        if self.transposed:
            txt = "[ ["
            for r in range(self.size[1]):
                for c in range(self.size[0]):
                    txt += "{{:{:s}}}".format(fmt).format(self[r,c]) # <- idx calculated for a trasposed matrix
                    if c < (self.size[0]-1):
                        txt += ", "
                if r < (self.size[1]-1):
                    txt += "], ["
                else:
                    txt += "]"
            txt += " ]"
        else:
            txt = "[ ["
            for r in range(self.size[0]):
                for c in range(self.size[1]):
                    txt += "{{:{:s}}}".format(fmt).format(self[r,c])
                    if c < (self.size[1]-1):
                        txt += ", "
                if r < (self.size[0]-1):
                    txt += "], ["
                else:
                    txt += "]"
            txt += " ]"

        return txt

    def __repr__(self):
        return "<class kMatrix {:s}>".format(self._do_format("f"))

    def __format__(self, fmt):
        return self._do_format(fmt)

    @property
    def T(self):
        return kMatrix(self.vector, size=self.size, transposed=True)

    #( --- iter --- )#
    def __iter__(self):
        if self.transposed:
            R = self.size[1]
            C = self.size[0]
        else:
            R = self.size[0]
            C = self.size[1]

        for r in range(R):
            for c in range(C):
                yield self[r,c]

    #( --- sum --- )#
    def __add__(self, y):
        if isinstance(y, kMatrix):
            assert self.size == y.size

            # the result shall have the right size, but not transposed
            # 1) A + B:
            if not self.transposed and not y.transposed:
                ret = kMatrix( [i+j for i,j in zip(self, y)], size=self.size)
            # 2) A + B.T:
            # 3) A.T + B:
            elif self.transposed and (not y.transposed):
                # self is flagged with "transposed", and the indexes will be calculated accordingly;
                # therefore, the 'tmp' matrix does not to be flagged again.
                print("self (c.T) 3x2")
                print(self)
                tmp = kMatrix(self)
                print("self.T")
                print(tmp)
                for i in tmp:
                    print("{:d} ".format(i))
                print("y")
                print(y)
                for i in y:
                    print("{:d} ".format(i))

                ret = kMatrix([i+j for i,j in zip( kMatrix(self, transposed=True), y )], size=y.size)
            # 4) A.T + B.T:
        elif isinstance(y, float) or isinstance(y, int):
            ret = kMatrix( [i+y for i in self.vector], size=self.size, transposed=self.transposed )
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

    def __eq__(self, y):
        tol = 1e-10
        ret = True
        if isinstance(y, kMatrix):
            if (self.transposed != y.transposed) or (self.size != y.size):
                ret = False

            if 1 == 1:
                pass

            if not all( [ abs(i-j) <= max( [abs(i), abs(j)] )*tol for i,j in zip(self.vector, y.vector) ] ):
                ret = False
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
if __name__ == "__main__":
    print("==== __init__() ====")
    a = kMatrix([i for i in range(6)], size=(3,2))

    print("==== __getitem__() ====")
    val = 0
    for r in range(3):
        for c in range(2):
            assert a[r,c] == val
            val += 1

    print("a = {:f}".format(a))
    print("a.T = {:f}".format(a.T))

    val = [0,2,4,1,3,5]
    for r in range(2):
        for c in range(3):
            assert a.T[r,c] == val.pop(0)

    print("==== iter ====")
    a = kMatrix([i for i in range(6)], size=(2,3)) # 0..5 [2x3]
    b = kMatrix(a, transposed=True) # 0..5 [3x2]

    # classic:
    val = 0
    for r in range(2):
        for c in range(3):
            print("a({:d},{:d}) = {:d}".format(r,c,a[r,c]))
            assert a[r,c] == val
            val += 1

    # iter():
    val = 0
    for i in a:
        assert i == val
        val += 1

    #classic:
    print()
    val = [0,3,1,4,2,5]
    for r in range(3):
        for c in range(2):
            print("b({:d},{:d}) = {:d}".format(r,c,b[r,c]))
            assert b[r,c] == val.pop(0)

    # iter() (transposed):
    val = [0,3,1,4,2,5]
    for i in b:
        assert i == val.pop(0)

    print("==== sum ====")
    a = kMatrix([i for i in range(6)], size=(3,2)) # 0..5
    b = kMatrix([i+1 for i in range(6)], size=(3,2)) # 1..6
    c = kMatrix(a, size=(2,3)) # 0..5

    assert a+b == kMatrix( [1,3,5,7,9,11], size=(3,2), transposed=False )
    assert (a+b).transposed == False

    assert c.T+a == kMatrix( [0,4,3,7,6,10], size=(3,2) )
    assert (c.T+a).transposed == False

    assert a+(-1) == kMatrix( [0,2,4,6,8,10], size=(3,2), transposed=False )
    assert a.T+b == None
    assert a.T+b.T == kMatrix( [0,3,5,7,9,11], size=(3,2), transposed=True )
    print("a + b    = {:f}".format(a+b))
    print("a + (-1) = {:f}".format(a+(-1)))
    print("(-1) + a = {:f}".format((-1)+a))
    a += 7
    assert a == kMatrix( [7,8,9,10,11,12], size=(3,2), transposed=False)
    print("a += 7: {:f}".format(a))
    print("a.T + b: {:f}".format(a.T+b))

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
