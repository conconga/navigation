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
#import numpy as np
#import copy
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kMatrixBuffer:
    def __init__(self, content):
        """
        The content shall be something convertable to a list, which will
        be stored.
        """

        self.buffer = list(content)
        self.len    = len(self.buffer)

    def __getitem__(self, index):
        assert index < self.len
        return self.buffer[index]

    def __setitem__(self, index, val):
        self.buffer[index] = val

    def __len__(self):
        return self.len

    def __str__(self):
        txt  = "buffer      = {:s}\n".format(self.buffer.__str__())
        txt += "len(buffer) = {:d}".format(self.len)
        return txt

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kMatrixInfo:
    def __init__(self, size, isTransposed):
        """
        The size of the original matrix, not transposed, is preserved and
        shall not be adjusted when the matrix is converted to a transposed one.
        """
        self.size_orig    = size
        self.isTransposed = isTransposed

    def __str__(self):
        txt  = "size_orig   = {:s}\n".format(self.size_orig.__str__())
        txt += "isTransp    = {:d}".format(self.isTransposed)
        return txt

    @property
    def size_rotated(self):
        return (self.size_orig[1], self.size_orig[0])

    @property
    def length(self):
        return self.size_orig[0] * self.size_orig[1]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kMatrixIdx:
    def __init__(self, info):
        """
        asTransposed: [True/False] defines how to calculate the indexing.
        """
        self.size_orig    = info.size_orig
        self.asTransposed = info.isTransposed

    def __getitem__(self, index_rc):
        return self._get_index(index_rc[0], index_rc[1])

    def _get_index(self, row, column):
        """
        For example:
        A matrix [2x3] is created with
        Buffer = [ a,b,c,d,e,f ]
        The matrix is then:
        [ [a,b,c],
          [d,e,f] ]
        The element 'e' at position (1,1) has index 4, and the element 'c'
        at position (0,2) has index 2.
        c(0,2): 0x3 + 2 = 2
        e(1,1): 1x3 + 1 = 4

        Now the buffer is used by a transposed matrix [3x2].
        [ [a,d],
          [b,e],
          [c,f] ]
        The element 'c' is now (2,0), but its index at the buffer remains 2.
        c(2,0): 2 + 0x2 = 2
        e(0,2): 0 + 2x2 = 4
        """

        if self.asTransposed:
            idx = row + (column * self.size_orig[1])
        else:
            idx = (row * self.size_orig[1]) + column

        return idx

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
class kMatrix:
    """
    The matrix is stored as a vector with the elements fill rows before columns.
    """
    def __init__(self, val, transposed=False, size=None):
        # loading:
        if isinstance(val, kMatrix):
            self.buffer = val.buffer # no deepcopy !!
            if size is None:
                size = val.info.size_orig

        elif isinstance(val, list) or isinstance(val, tuple):
            self.buffer = kMatrixBuffer(val)
            if size is None:
                raise(NameError("the size must be specified here"))

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        self.info = kMatrixInfo(size, transposed)
        self.idx  = kMatrixIdx(self.info)
        assert len(self.buffer) == self.info.length

    def __getitem__(self, index):
        """
        index = (row, column)
        """
        return self.buffer[self.idx[index]]

    def __setitem__(self, index, val):
        """
        index = (row, column)
        """
        self.buffer[self.idx[index]] = val

    def print_summary(self):
        print("vvv sumary vvv")
        print(self.info)
        print(self.buffer)
        print("^^ sumary ^^")

    def _get_RC(self):
        if self.info.isTransposed:
            R,C = self.info.size_rotated
        else:
            R,C = self.info.size_orig
        return R,C

    def _do_format(self, fmt):
        #self.print_summary()

        R,C = self._get_RC()

        txt = "[ ["
        for r in range(R):
            for c in range(C):
                txt += "{{:{:s}}}".format(fmt).format(self[r,c])
                if c < (C-1):
                    txt += ", "
            if r < (R-1):
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
        if self.info.isTransposed:
            return kMatrix(self, size=self.info.size_orig, transposed=False)
        else:
            return kMatrix(self, size=self.info.size_orig, transposed=True)

    #( --- iter --- )#
    def __iter__(self):
        R,C = self._get_RC()
        for r in range(R):
            for c in range(C):
                yield self[r,c]

    #( --- sum --- )#
    def __add__(self, y):
        if isinstance(y, kMatrix):
            # the result shall have the right size, but not transposed

            # 1) A + B:
            if (not self.info.isTransposed) and (not y.info.isTransposed):
                assert self.info.size_orig == y.info.size_orig
                ret = kMatrix( [i+j for i,j in zip(self, y)], size=self.info.size_orig)

            # 2) A + B.T:
            elif (not self.info.isTransposed) and (y.info.isTransposed):
                assert self.info.size_orig == y.info.size_rotated
                ret = kMatrix([i+j for i,j in zip( self, kMatrix(y, transposed=True) )], size=self.info.size_orig )

            # 3) A.T + B:
            elif self.info.isTransposed and (not y.info.isTransposed):
                assert self.info.size_rotated == y.info.size_orig
                ret = kMatrix([i+j for i,j in zip( kMatrix(self, transposed=True), y )], size=y.info.size_orig )

            # 4) A.T + B.T:
            elif (self.info.isTransposed) and (y.info.isTransposed):
                assert self.info.size_rotated == y.info.size_rotated
                ret = kMatrix([i+j for i,j in zip( kMatrix(self, transposed=True), kMatrix(y, transposed=True) )],
                              size=self.info.size_rotated)

        elif isinstance(y, float) or isinstance(y, int):
            ret = kMatrix( [i+y for i in self], size=self.info.size_orig, transposed=self.info.isTransposed )
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))
        return ret

    def __radd__(self, y):
        return self.__add__(y)

    def __iadd__(self, y): # +=
        q = self.__add__(y)
        self.buffer = q.buffer
        return(self)

    #( --- negative signal --- )#
    def __neg__(self):
        return kMatrix( [-i for i in self], size=self.info.size_orig )

    #( --- difference --- )#
    def __sub__(self, y):
        if isinstance(y, kMatrix):
            # the result shall have the right size, but not transposed

            # 1) A - B:
            if (not self.info.isTransposed) and (not y.info.isTransposed):
                assert self.info.size_orig == y.info.size_orig
                ret = kMatrix( [i-j for i,j in zip(self, y)], size=self.info.size_orig)

            # 2) A - B.T:
            elif (not self.info.isTransposed) and (y.info.isTransposed):
                assert self.info.size_orig == y.info.size_rotated
                ret = kMatrix([i-j for i,j in zip( self, kMatrix(y, transposed=True) )], size=self.info.size_orig )

            # 3) A.T - B:
            elif self.info.isTransposed and (not y.info.isTransposed):
                assert self.info.size_rotated == y.info.size_orig
                ret = kMatrix([i-j for i,j in zip( kMatrix(self, transposed=True), y )], size=y.info.size_orig )

            # 4) A.T - B.T:
            elif (self.info.isTransposed) and (y.info.isTransposed):
                assert self.info.size_rotated == y.info.size_rotated
                ret = kMatrix([i-j for i,j in zip( kMatrix(self, transposed=True), kMatrix(y, transposed=True) )],
                              size=self.info.size_rotated)

        elif isinstance(y, float) or isinstance(y, int):
            ret = kMatrix( [i-y for i in self], size=self.info.size_orig, transposed=self.info.isTransposed )
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))
        return ret

    def __rsub__(self, y): # y - self
        return y + (-self)

    def __isub__(self, y): # -=
        q = self.__sub__(y)
        self.buffer = q.buffer
        return(self)

    #( --- multiplication --- )#
    def __mul__(self, y):
        if isinstance(y, kMatrix):
            Rs,Cs = self._get_RC()
            Ry,Cy = y._get_RC()
            assert Cs == Ry
            ret = kMatrix( [0.0 for i in range(Rs*Cy)], size=(Rs, Cy), transposed=False )

            for i in range(Rs):
                for j in range(Cy):
                    val = 0.0
                    for k in range(Cs):
                        val += self[i,k] * y[k,j]
                    ret[i,j] = val

        elif isinstance(y, int) or isinstance(y, float):
            ret = kMatrix( [y*i for i in self], size=self.info.size_orig, transposed=False )

        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

    def __rmul__(self, y):
        return self.__mul__(y)

    def __imul__(self, y): # *=
        q = self.__mul__(y)
        self.buffer = q.buffer
        return(self)

    #( --- division --- )#
    def __truediv__(self, y):
        raise(NameError("WTF??"))

    def __rtruediv__(self, y):
        raise(NameError("WTF??"))

    def __itruediv__(self, y): # /=
        raise(NameError("WTF??"))

    #( --- miscelaneous --- )#
    def __eq__(self, y):
        tol = 1e-10
        ret = True
        if isinstance(y, kMatrix):
            if (self.info.isTransposed != y.info.isTransposed) or (self.info.size_orig != y.info.size_orig):
                ret = False

            if 1 == 1:
                pass

            if not all( [ abs(i-j) <= max( [abs(i), abs(j)] )*tol for i,j in zip(self, y) ] ):
                ret = False
        else:
            raise(NameError("not prepared for type '{:s}'".format(str(type(y)))))

        return ret

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
if __name__ == "__main__":
    print("==== __init__() ====")
    a = kMatrix([i for i in range(6)], size=(3,2))

    print("==== no deepcopy ====")
    b = kMatrix(a, size=(3,2))
    a.buffer.buffer[0] = -10
    assert b[0,0] == -10

    print("==== __getitem__() ====")
    a = kMatrix([i for i in range(6)], size=(3,2))
    val = 0
    for r in range(3):
        for c in range(2):
            assert a[r,c] == val
            val += 1

    print("a = {:f}".format(a))
    print("a.T = {:f}".format(a.T))
    print("a.T.T = {:f}".format(a.T.T))
    print("a.T.T.T = {:f}".format(a.T.T.T))

    val = [0,2,4,1,3,5]
    for r in range(2):
        for c in range(3):
            assert a.T[r,c] == val.pop(0)

    print("==== transpose ====")
    assert a.T.T == a
    assert a.T.T.T == a.T

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
    assert c.T+a == kMatrix( [0,4,3,7,6,10], size=(3,2), transposed=False )
    assert b+c.T == kMatrix( [1,5,4,8,7,11], size=(3,2), transposed=False )
    assert a.T + b.T == kMatrix( [1,5,9,3,7,11], size=(2,3), transposed=False )

    assert a+(-1) == kMatrix( [-1,0,1,2,3,4], size=(3,2), transposed=False )
    assert (-1)+a == kMatrix( [-1,0,1,2,3,4], size=(3,2), transposed=False )

    a += 7
    assert a == kMatrix( [7,8,9,10,11,12], size=(3,2), transposed=False)
    print("a += 7: {:f}".format(a))

    print("==== negative signal ====")
    assert -a == kMatrix( [-7,-8,-9,-10,-11,-12], size=(3,2) )
    print("-a = {:f}".format(-a))

    print("==== difference ====")
    a = kMatrix([i for i in range(6)], size=(3,2)) # 0..5
    b = kMatrix([i+1 for i in range(6)], size=(3,2)) # 1..6
    c = kMatrix(a, size=(2,3)) # 0..5

    assert a-b      == kMatrix( [-1 for i in range(6) ], size=(3,2), transposed=False )
    assert a.T-c    == kMatrix( [0,1,2,-2,-1,0], size=(2,3), transposed=False )
    assert a-c.T    == kMatrix( [0,-2,1,-1,2,0], size=(3,2), transposed=False )
    assert a.T-b.T  == kMatrix( [-1 for i in range(6)], size=(2,3), transposed=False )

    assert a-1 == kMatrix( [-1,0,1,2,3,4], size=(3,2), transposed=False )
    assert 1-a == kMatrix( [1,0,-1,-2,-3,-4], size=(3,2), transposed=False )

    a -= 3.0
    assert a == kMatrix( [-3,-2,-1,0,1,2], size=(3,2), transposed=False )
    print("a -= 3.0: {:f}".format(a))

    print("==== multiplication ====")
    a = kMatrix([i for i in range(6)], size=(3,2)) # 0..5
    b = kMatrix([i+1 for i in range(6)], size=(3,2)) # 1..6
    c = kMatrix(a, size=(2,3)) # 0..5

    #import numpy as np
    #A = np.asarray([[0,1],[2,3],[4,5]])
    #B = np.asarray([[1,2],[3,4],[5,6]])
    #C = np.asarray([[0,1,2],[3,4,5]])
    #print(A.dot(B.T))
    #print(A.T.dot(B))
    #print(A.T.dot(C.T))
    #print(A.dot(C))

    assert a*b.T == kMatrix([2,4,6,8,18,28,14,32,50], size=(3,3), transposed=False)
    assert a.T*b == kMatrix([26,32,35,44], size=(2,2), transposed=False)
    assert a.T*c.T == kMatrix([10,28,13,40], size=(2,2), transposed=False)
    assert a*c == kMatrix([3,4,5,9,14,19,15,24,33], size=(3,3), transposed=False)

    assert a*2.0 == kMatrix( [0,2,4,6,8,10], size=(3,2), transposed=False)
    assert 2.0*b == kMatrix( [2,4,6,8,10,12], size=(3,2), transposed=False)
    b *= -1.0
    assert b == kMatrix( [-1,-2,-3,-4,-5,-6], size=(3,2), transposed=False)

    try:
        print(a*a)
    except:
        print("ok")
        a = None

    if a is not None:
        raise(NameError("it should not reach here"))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
