#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: kArray.py
Beschreibung: Mathematical manipulation of vectors and matrices.
Autor: Luciano Auguto Kruk
Erstellt am: 06.09.2025
Version: 1.0.0
Lizenz: 
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
import numpy as np
import math
#import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArrayCommon:

    TYPE_ARRAY       = 0
    TYPE_VERTICAL    = 1
    TYPE_HORIZONTAL  = 2
    TYPE_SINGLEVALUE = 3

    def _type(self, val):
        """
            TYPE_ARRAY       = 2D array
            TYPE_VERTICAL    = vertical vector
            TYPE_HORIZONTAL  = horizontal vector
            TYPE_SINGLEVALUE = single value array
        """

        assert isinstance(val, np.ndarray)
        assert all( [i>0 for i in val.shape] )

        size = val.shape
        if len(size) == 1:
            if size[0] == 1:
                return self.TYPE_SINGLEVALUE
            else:
                return self.TYPE_HORIZONTAL
        elif len(size) == 2:
            if size == (1,1):
                return self.TYPE_SINGLEVALUE
            if size[0] == 1:
                return self.TYPE_HORIZONTAL
            elif size[1] == 1:
                return self.TYPE_VERTICAL
            else:
                return self.TYPE_ARRAY
        else:
            print("::error::")
            print(val)
            raise(NameError("I am not prepared for this array"))


    #def _do_format_1D(self, fmt, C):
    #    txt = "[ "
    #    for c in range(C):
    #        txt += "{{:{:s}}}".format(fmt).format(self.array[c])
    #        if c < (C-1):
    #            txt += ", "
    #    txt += " ]{:s}".format("" if self.hvector else ".T")
    #    return txt

    def _do_format_2D(self, fmt, R, C):
        txt = "[ ["
        for r in range(R):
            for c in range(C):
                txt += "{{:{:s}}}".format(fmt).format(self.array[r,c])
                if c < (C-1):
                    txt += ", "
            if r < (R-1):
                txt += "], ["
            else:
                txt += "]"
        txt += " ]"
        return txt

    def _do_format(self, fmt):
        size  = self.array.shape
        vtype = self._type(self.array)
        txt = self._do_format_2D(fmt, *size)
        return txt

    #( --- indexing ---)#
    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.array[idx]
        else:
            return self.array[*idx]

    def __setitem__(self, idx, val):
        if isinstance(idx, int):
            self.array[idx] = val
        else:
            self.array[*idx] = val

    #( --- miscelaneous --- )#
    def __eq__(self, y):
        tol = 1e-10
        ret = True
        y   = kArray(y)

        if self.array.shape != y.array.shape:
            ret = False

        if not all( [ abs(i-j) <= max( [abs(i), abs(j)] )*tol for i,j in zip(self, y) ] ):
            ret = False

        return ret

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kArray (kArrayCommon):
    def __init__(self, val, hvector=None):
        """
        When 'val' is given as a list, 'hvector' is used to indicate whether the
        vector to be created is horizontal (row, True) or vertical (column, False).

        if 'hvector' is None, then the selection is automatic when possible.
        """

        if isinstance(val, (list, tuple)):
            val = np.asarray( val )
        elif isinstance(val, (int, float)):
            val = np.asarray( [val] )
        elif isinstance(val, kArray):
            val = val.array

        shape = val.shape
        assert 1 <= len(shape) <= 2
        vtype = self._type(val)

        if vtype == self.TYPE_ARRAY:
            self.array = val

        elif vtype in [ self.TYPE_HORIZONTAL, self.TYPE_VERTICAL ]:
            if hvector == True:
                self.array = val.squeeze().reshape(1,-1)
            elif hvector == False:
                self.array = val.squeeze().reshape(-1,1)
            else: # hvector==None
                if vtype == self.TYPE_HORIZONTAL:
                    self.array = val.squeeze().reshape(1,-1)
                else:
                    self.array = val.squeeze().reshape(-1,1)

        elif vtype == self.TYPE_SINGLEVALUE:
            self.array = val.squeeze().reshape(1,1)

        else:
            print("::error::")
            print(val)
            print(vtype)
            raise(NameError("what is this?"))

        #print("self.array =")
        #print(self.array)

    def __repr__(self):
        #txt = "<class {:s} ".format(str(self.__class__))
        #return txt + "{:s}>".format(self._do_format("f"))
        return "{:s} = !{:s}!".format(str(self.__class__), self._do_format("f"))


    def __format__(self, fmt):
        return self._do_format(fmt)

    @property
    def T(self):
        vtype = self._type(self.array)
        if vtype == self.TYPE_HORIZONTAL:
            return self.__class__( self.array, hvector=False )
        elif vtype == self.TYPE_VERTICAL:
            return self.__class__( self.array, hvector=True )
        else:
            return self.__class__( self.array.T )

    #( --- iter --- )#
    def __iter__(self):
        if self._type(self.array) == self.TYPE_SINGLEVALUE:
            yield float(self.array.squeeze())
        else:
            temp = self.array.reshape(1,-1).squeeze()
            for i in temp:
                yield i

    #( --- sum --- )#
    def __add__(self, y):
        if isinstance(y, kArray):
            if self.array.shape != y.array.shape:
                raise(NameError("both arrays shall have the same dimensions"))
            else:
                ret = self.__class__( self.array + y.array )

        elif isinstance(y, int) or isinstance(y, float):
            ret = self.__class__( y + self.array )

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

    def __radd__(self, y): # y + self
        return self.__add__(y)

    def __iadd__(self, y): # +=
        self.array += y
        return(self)

    #( --- negative signal --- )#
    def __neg__(self):
        return self.__class__( -self.array )

    #( --- difference --- )#
    def __sub__(self, y):
        if isinstance(y, kArray):
            if self.array.shape != y.array.shape:
                raise(NameError("both vector shall have the same dimensions"))
            else:
                ret = self.__class__( self.array - y.array )

        elif isinstance(y, int) or isinstance(y, float):
            ret = self.__class__( self.array - y )

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

    def __rsub__(self, y):
        # y - self
        return (-self).__add__(y)

    def __isub__(self, y): # -=
        self.array = self.array - y
        return(self)

    #( --- multiplication --- )#
    def __mul__(self, y):
        if isinstance(y, kArray):
            axb   = self.array.dot(y.array)
            vtype = self._type( np.asarray(axb) )
            if vtype == self.TYPE_HORIZONTAL:
                ret = self.__class__( axb, hvector=True )
            elif vtype == self.TYPE_VERTICAL:
                ret = self.__class__( axb, hvector=False )
            else:
                ret = self.__class__( axb )

        elif isinstance(y, int) or isinstance(y, float):
            ret = self.__class__( self.array * y )

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

    def __rmul__(self, y):
        # y * self
        return self.__mul__(y)

    def __imul__(self, y): # *=
        self = self * y
        return(self)

    #( --- miscelaneous --- )#
    def __abs__(self):
        vtype = self._type(self.array)
        if vtype in [ self.TYPE_HORIZONTAL, self.TYPE_VERTICAL ]:
            return math.sqrt( sum( [i**2 for i in self.array.squeeze()] ))
        elif vtype == self.TYPE_SINGLEVALUE:
            return abs(self.array.squeeze())
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

    def to_list(self):
        return self.array.reshape(1,-1).squeeze().tolist()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def tests_general():
    print("==== type() ====")

    a = kArray([[1,2],[3,4]]) # only to get in with a valid object
    assert a._type( np.asarray( [1] )) == a.TYPE_SINGLEVALUE
    assert a._type( np.asarray( [1,2,3] )) == a.TYPE_HORIZONTAL
    assert a._type( np.asarray( [[1,2,3]] )) == a.TYPE_HORIZONTAL
    assert a._type( np.asarray( [[1]] )) == a.TYPE_SINGLEVALUE
    assert a._type( np.asarray( [[1,2]] )) == a.TYPE_HORIZONTAL
    assert a._type( np.asarray( [[1],[2]] )) == a.TYPE_VERTICAL
    assert a._type( np.asarray( [[1,2],[3,4]] )) == a.TYPE_ARRAY

    try:
        ok = False
        a._type( np.asarray( [] ))
    except:
        ok = True
    if not ok:
        raise(NameError("Error"))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def tests_vector():
    print("==== __init__() ====")

    print("kArray([1]) = ")
    print(kArray([1]))
    print("kArray([1,2,3]) = ")
    print(kArray([1,2,3]))
    print("kArray([1,2,3],True) = ")
    print(kArray([1,2,3],True))
    print("kArray([1,2,3],False) = ")
    print(kArray([1,2,3],False))
    print("kArray([[1,2,3]], False) = ")
    print(kArray([[1,2,3]], False))
    print("kArray([[1],[2]], True) = ")
    print(kArray([[1],[2]], True))

    print("==== transpose ====")
    a = kArray([1,2,3], hvector=True)
    print("a = {:f}".format(a))
    print("a.T = {:f}".format(a.T))
    print("a.T.T = {:f}".format(a.T.T))
    print("a.T.T.T = {:f}".format(a.T.T.T))
    print("a = {:f}".format(a))

    print("==== iter ====")
    a = kArray( [1,2,3], hvector=False )

    # iter():
    for i,j in zip(a,[1,2,3]):
        assert i == j

    print("==== eq ====")
    a = kArray([1,2,3])
    b = kArray([2,2,3])
    assert a == a
    assert a != b
    print("a == a: {:s}".format((a==a).__str__()))
    print("a == b: {:s}".format((a==b).__str__()))

    print("==== sum ====")
    a = kArray((1,2,3))
    b = kArray((4,5,6))
    assert a+b == kArray([5,7,9])
    assert a+(-1) == kArray([0,1,2])
    assert (-1)+a == kArray([0,1,2])
    print("a + b    = {:f}".format(a+b))
    print("a + (-1) = {:f}".format(a+(-1)))
    print("(-1) + a = {:f}".format((-1)+a))
    a += 7
    assert a == kArray([8,9,10])
    print("a += 7: {:f}".format(a))

    print("==== negative signal ====")
    print("-a = {:f}".format(-a))
    assert -a == kArray([-8,-9,-10])

    print("==== difference ====")
    a = kArray((1,2,3))
    b = kArray((4,5,6))
    assert a-b == kArray([-3, -3, -3])
    assert 3.-a == kArray([2,1,0])
    print("a - b   = {:f}".format(a-b))
    print("3.0 - a = {:f}".format(3.0-a))
    a -= 3.0
    assert a == kArray([-2,-1,0])
    print("a -= 3.0: {:f}".format(a))

    print("==== multiplication ====")
    a = kArray([1,2,3])
    b = kArray((4,5,6))
    try:
        print("a * b = {:f}".format(a*b))
    except:
        print("a * b =")
        print("    ^^^ERROR")

    assert a*b.T == 32
    assert 7.*a == kArray([7,14,21])
    print("a.T     = {:f}".format(a.T))
    print("7.0 * a = {:f}".format( 7.0 * a ))
    a *= 10
    assert a == kArray([10,20,30])
    print("a *= 10 : {:f}".format(a))

    print("==== norm ====")
    a = kArray([1,1])
    assert abs(abs(a) - math.sqrt(2)) < 1e-10
    print("||a|| = {:f}".format(abs(a)))
    a = kArray(-32)
    assert abs(abs(a) - 32.0) < 1e-10

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
def tests_matrix():
    print("==== __init__() ====")
    print("kArray( [[1,2,3],[4,5,6]] ) =")
    print(kArray( [[1,2,3],[4,5,6]] ))

    print("==== transpose ====")
    a = kArray( [[1,2],[3,4]] )
    print("a = {:f}".format(a))
    print("a.T = {:f}".format(a.T))
    print("a.T.T = {:f}".format(a.T.T))
    print("a.T.T.T = {:f}".format(a.T.T.T))

    val = [1,3,2,4]
    for i in a.T:
        assert i == val.pop(0)

    assert a.T.T == a
    assert a.T.T.T == a.T

    print("==== sum ====")
    a = kArray( [[1,2,3], [4,5,6]] )
    b = kArray( [[0,1,2], [3,4,5]] )
    c = kArray( [[10,11],[12,13],[14,15]] )

    assert a+b == kArray( [[1,3,5], [7,9,11]] )
    assert c.T+a == kArray( [[11,14,17], [15,18,21]] )
    assert b+c.T == kArray( [[10,13,16], [14,17,20]] )

    assert a+(-1) == b
    assert (-1)+a == b

    a += 7
    assert a == kArray( [[8,9,10], [11,12,13]] )
    print("a += 7: {:f}".format(a))

    print("==== negative signal ====")
    a = kArray( [[1,2,3], [4,5,6]] )
    assert -a == kArray( [[-1, -2, -3], [-4, -5, -6]] )
    print("-a = {:f}".format(-a))

    print("==== difference ====")
    a = kArray( [[1,2,3], [4,5,6]] )
    b = kArray( [[0,1,2], [3,4,5]] )
    c = kArray( [[10,11],[12,13],[14,15]] )

    assert a-b      == kArray( 1+np.zeros((2,3)) )
    assert a.T-c    == kArray( [[-9,-7],[-10,-8],[-11,-9]] )
    assert a-c.T    == kArray( [[-9,-10,-11],[-7,-8,-9]] )
    assert a.T-b.T  == kArray( np.ones((3,2)) )

    assert a-1 == kArray( [[0,1,2],[3,4,5]] )
    assert 1-a == kArray( [[0,-1,-2],[-3,-4,-5]] )

    a -= 3.0
    assert a == kArray( [[-2,-1,0],[1,2,3]] )
    print("a -= 3.0: {:f}".format(a))

    print("==== multiplication ====")
    a = kArray( [[1,2,3], [4,5,6]] ) # 2x3
    b = kArray( [[0,1,2], [3,4,5]] ) # 2x3
    c = kArray( [[10,11],[12,13],[14,15]] ) # 3x2
    d = kArray( [1,2,3], hvector=False ) # 3x1
    e = kArray( [2,3], hvector=True ) # 1x2

    assert a*b.T == kArray( [[8,26],[17,62]] )
    assert a.T*b == kArray( [[12,17,22],[15,22,29],[18,27,36]] )
    assert a*c == kArray( [[76, 82],[184,199]] )
    assert a*d == kArray( [[14],[32]] )
    assert e*a*d == 124.

    assert a*2.0 == kArray( [[2,4,6],[8,10,12]] )
    assert 2.0*b == kArray( [[0,2,4],[6,8,10]] )
    b *= -1.0
    assert b == kArray( [[0,-1,-2],[-3,-4,-5]] )

    try:
        ok = False
        print(a*a)
    except:
        ok = True

    if not ok:
        raise(NameError("it should not reach here"))

    print("==== indexing ====")
    a = kArray( [[1,2,3], [4,5,6]] ) # 2x3
    assert a[1,1] == 5
    assert list(a[1]) == [4,5,6]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

if __name__ == "__main__":
    tests_general()
    tests_vector()
    tests_matrix()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
