import math, os, gc, copy
import numpy as np
from ctypes import *
from numpy import pi

from ..tools.check import ispositiveint
from ..general import FQHE, abpath, mkdir
from ..core import interaction
from ..core.Hamiltonian import csr, csr_fqhe
from ..core.eigenvector import eigenvector as eigvec
#from ..core.eigenvector import eigenvector_k1
from ..core.Hilbertspace.SLFermion import SLFermion_wrapper as cfl

class FermionSL(FQHE):
    def __init__(self, n, n_o, num_thread=8):
        super().__init__(n, n_o, num_thread)
        self.__level = 0

        self.__A_2body = 0
        self.__A_3body = 0
        cfl.init_FQHE(self.n_o, self.n)

    def __repr__(self):
        return str(self.geometry)

    def clA(self):
        self.__A_2body = 0
        self.__A_3body = 0

    @property
    def geometry(self):
        gd = {}
        gd['n'] = self.n
        gd['n_phi'] = (self.n_o-1)/2
        gd['n_o'] = self.n_o
        gd['level'] = self.__level
        gd['num_thread'] = self.num_thread
        return gd

    @property
    def inter(self):
        return self.__inter

    @property
    def Lz_dim(self):
        return cfl.get_Lz_dim()

    def setinteraction(self, inter, calQ=False):
        if not isinstance(inter, interaction.interaction):
            raise ValueError("Variable inter is not an interaction instance!")
        self.__inter = inter
        self.__twobody = inter.twobody
        self.__A_2body = 0
        if not calQ:
            return
        geom = self.geometry
        #for ele in self.__twobody:
            #self.__A_2body += getattr(interaction, 'twobody'+ele['type'])(ele, geom)
        self.__A_2body = np.zeros((self.n_o,self.n_o,self.n_o,self.n_o), dtype=np.float64)
        cfl.createA2BdoywithNo(self.__twobody, self.__A_2body, self.n_o, self.num_thread)
        return self.__A_2body

    def copy(self):
        return copy.deepcopy(self)

    def Lz_space(self, Lz=0, copy=False):
        #if 0 > Lz :
        #    raise ValueError("k1 = ", k1,"N_phi = ", self.n_phi)
        cfl.create_Lz_space( round(Lz+self.n*(self.n_o-1)//2) )
        self.__Lz = Lz
        if copy:
            basis = np.zeros((cfl.get_Lz_dim(), self.n), dtype=np.int8)
            cfl.copy_basis(basis)
            return basis

    def Hilbertspace(self, Lz=0):
        #if 0 > k1 or k1 >= self.n_phi:
        #    raise ValueError("k1 = ", k1,"N_phi = ", self.n_phi)
        if (self.n_o!=cfl.n_o()) or (self.n!=cfl.n_e()):
            print("init!")
            print((self.n_o,cfl.n_o()),(self.n,cfl.n_e()))
            cfl.init_FQHE(self.n_o, self.n)
        self.Lz_space(Lz)

    def twobodyH(self, Lz=0, HermitianQ=0):
        if isinstance(self.__A_2body, int):
            self.setinteraction(self.__inter, calQ=True)
        self.Hilbertspace(Lz)
        cfl.setHermitian(HermitianQ)
        indptr = np.zeros(cfl.get_Lz_dim()+1,dtype=np.uint32, order = 'C')
        nnz = cfl.twoBodyCount(indptr, self.num_thread)
        gc.collect()
        indices = np.zeros(nnz, dtype=np.uint32, order = 'C')
        data = np.zeros(nnz, dtype=np.complex128, order = 'C')
        cfl.twoBodyGenerate(self.__A_2body, indptr, indices, data, self.num_thread)
        return csr_fqhe(indptr, indices, data, HermitianQ, self, Lz, "twobody", self.__inter)

class fqhe_sl_solver():
    '''
        instance(ins, FQHE)
        instance(inter, interaction)
        ["twobody"]
    '''
    def __init__(self, ins, inter, itype, Lz=0):
        if not isinstance(ins, FermionSL):
            raise ValueError(ins, "is not a FermionSL isinstance!")
        self.__Lz = Lz
        self.__fqhe = ins
        if itype not in ["twobody"]:
            raise ValueError(itype)
        self.__itype = itype
        self.__twobody = inter.twobody.copy()
        #self.__threebody = inter.threebody.copy()
        self.__fqhe.setinteraction(inter, calQ=False)

    def set_Lz(self, Lz):
        self.__Lz = Lz

    def setinteraction(self, inter, calQ=False):
        if not isinstance(inter, interaction.interaction):
            raise ValueError("Variable inter is not an interaction instance!")
        self.__inter = inter
        self.__twobody = inter.twobody
        #self.__threebody = inter.threebody
        cal = calQ
        self.__fqhe.setinteraction(inter, calQ=cal)

    @property
    def fqhe(self):
        return self.__fqhe

    @property
    def Lz(self):
        return self.__Lz

    @property
    def itype(self):
        return self.__itype

    @property
    def geometry(self):
        return self.__fqhe.geometry

    def eig(self, tol=10**(-15), vec=False, HermitianQ=0):
        csr = getattr(self.__fqhe, self.__itype+"H")(Lz=self.__Lz, HermitianQ=HermitianQ)
        ge, gv = csr.lanczos(tol = tol)
        if vec:
            return eigvec(gv, ge, self.__fqhe, self.k1, self.k2)
        else:
            return ge

    def eigs(self, tol=10**(-15), HermitianQ=0):
        csr = getattr(self.__fqhe, self.__itype+"H")(Lz=self.__Lz, HermitianQ=HermitianQ)
        ge = csr.lanczosIRL(tol=tol,k=4)[0]
        return ge

    def eigsv(self, k=4, tol=10**(-15), vec=True, HermitianQ=1):
        csr = getattr(self.__fqhe, self.__itype+"H")(Lz=self.__Lz, HermitianQ=HermitianQ)
        ge = csr.lanczosIRL(tol=tol,k=k)
        out = []
        if vec:
            for i in range(k):
                out.append(eigvec(np.transpose(ge[1])[i], ge[0][i], self.__fqhe, self.Lz))
            return out
        else:
            return ge[0]

    def Hamiltonian(self, HermitianQ=0):
        csr = getattr(self.__fqhe, self.__itype+"H")(Lz=self.__Lz, HermitianQ=HermitianQ)
        return csr



