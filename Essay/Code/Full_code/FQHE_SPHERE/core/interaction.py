import numpy as np
from scipy import sparse
from numpy import pi
from ctypes import *
import math, cmath, gc, copy, os
from itertools import product
from scipy.special import genlaguerre

from ..general import FQHE
#fermion = CDLL(os.path.dirname(__file__) + "/Hilbertspace/SLFermion/fermionSL.so")
#boson = CDLL(os.path.dirname(__file__) + "/threebodyinteraction/boson.so")
#void A_Three_Body_0(complex double *A_list。 complex double *A_list_0, double *qx, double *qy, \
#					int bound, int N_phi, int num_thread)
#boson.A_Three_Body_0.argtypes = [
#							np.ctypeslib.ndpointer(dtype=np.complex128, ndim=4,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.complex128, ndim=4,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.float64, ndim=2,flags='C_CONTIGUOUS'),
#							c_int, c_int, c_int, c_double, c_double, c_double]
#void A_Three_Body_0(complex double *A_list, complex double *A_list_0, double *qx, double *qy, \
#					int bound, int N_phi, int num_thread, double g11, double g12, double g22);
#fermion.A_Three_Body_0.argtypes = [
#							np.ctypeslib.ndpointer(dtype=np.complex128, ndim=4,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.complex128, ndim=4,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags='C_CONTIGUOUS'),
#							np.ctypeslib.ndpointer(dtype=np.float64, ndim=2,flags='C_CONTIGUOUS'),
#							c_int, c_int, c_int, c_double, c_double, c_double]


up_bound = 113

class interaction():
	'''
		#setcoulomb(self, inten=1.0, ani_inter=[0.0, 0.0], ani_mass=[0.0, 0.0])
		setpseudo2body(self, inten=1.0, order=-1, ani_intrin=[0.0, 0.0])
		#setpseudo3body(self, inten=1.0, anti=0.0)
	'''
	def __init__(self):
		self.__twobody = np.zeros(up_bound, dtype=np.float64)
		#self.__threebody = []

	#def setcoulomb(self, inten=1.0, ani_inter=[0.0, 0.0], ani_mass=[0.0, 0.0], level=0):
	#	if inten >= 0:
	#		self.__twobody.append({'type':'coulomb', 'inten':inten,
	#							 'ani_inter':ani_inter, 'ani_mass':ani_mass, 'level':level})
	#	else:
	#		raise ValueError("inten = ", inten, " < 0")

	#def setcoulomb_sqrt(self, inten=1.0, ani_inter=[0.0, 0.0], ani_mass=[0.0, 0.0], level=0, d=4):
	#	if inten >= 0:
	#		self.__twobody.append({'type':'coulomb_sqrt', 'inten':inten,
	#							 'ani_inter':ani_inter, 'ani_mass':ani_mass, 'level':level, 'd':d})
	#	else:
	#		raise ValueError("inten = ", inten, " < 0")

	#def setcoulomb_exp(self, inten=1.0, ani_inter=[0.0, 0.0], ani_mass=[0.0, 0.0], level=0, d=4):
	#	if inten >= 0:
	#		self.__twobody.append({'type':'coulomb_exp', 'inten':inten,
	#							 'ani_inter':ani_inter, 'ani_mass':ani_mass, 'level':level, 'd':d})
	#	else:
	#		raise ValueError("inten = ", inten, " < 0")

	def setpseudo2body(self, inten=[0.0, 1.0]):
		self.__twobody = np.zeros(up_bound, dtype=np.float64)
		for i in range(len(inten)):
			self.__twobody[i] = inten[i]

	#def setgpseudo2body(self, inten=1.0, orderm=-1, ordern=-1, ani_intrin=[0.0, 0.0]):
	#	if inten >= 0 and isinstance(orderm, int) and isinstance(ordern, int):
	#		self.__twobody.append({'type':'gpseudo2body', 'inten':inten,           
	#							 'orderm':orderm, 'ordern':ordern, 'ani_intrin':ani_intrin})
	#	else:
	#		raise ValueError("inten = ", inten, " <= 0")

	#def setpseudo3body(self, inten=1.0, anti=0.0, ani_intrin=[0.0, 0.0]):
	#	self.__threebody.append({'type':'threebody', 'inten':inten,
	#								'anti':anti, 'ani_intrin':ani_intrin})

	@property
	def twobody(self):
		return self.__twobody

	#@property
	#def threebody(self):
	#	return self.__threebody

def _twobodyq2(geom, ani):
	Q, phi = ani
	coord = np.arange(1 - up_bound, up_bound)
	qx = coord.reshape(2*up_bound-1, 1)/geom['lx']
	qy = coord/geom['ly'] - coord.reshape(2*up_bound-1, 1)*math.cos(geom['theta'])/geom['lx']
	qx = 2*pi*qx
	qy = 2*pi*qy/math.sin(geom['theta'])
	g11= np.cosh(Q)+np.cos(phi)*np.sinh(Q)
	g12= np.sin(phi)*np.sinh(Q)
	g22= np.cosh(Q)-np.cos(phi)*np.sinh(Q) #P.R.B 155140 (2018) Eq.(2)
	return g11*(qx**2)+g22*(qy**2)+2*g12*qx*qy

def _twobodyvq(geom):
	coord = np.arange(1 - up_bound, up_bound)
	qx = coord.reshape(2*up_bound-1, 1)/geom['lx']
	qy = coord/geom['ly'] - coord.reshape(2*up_bound-1, 1)*math.cos(geom['theta'])/geom['lx']
	qx = 2*pi*qx
	qy = 2*pi*qy/math.sin(geom['theta'])
	return (qx+1j*qy)/np.sqrt(2.0)

def twobodycoulomb(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")
	coord = np.arange(1 - up_bound, up_bound)

	q2inter = _twobodyq2(geom, inter['ani_inter'])
	q2mass = _twobodyq2(geom, inter['ani_mass'])
	q2inter[up_bound-1, up_bound-1] = pi
	temp = np.exp(-q2mass/2)
	temp = temp*(2*pi/np.sqrt(q2inter))/(2*geom['area'])
	temp[up_bound-1, up_bound-1] = 0
	aa = np.zeros(inter['level']+1)
	aa[inter['level']] = 1.0
	_interaction = np.polynomial.laguerre.Laguerre(aa)
	temp *= (_interaction(q2mass/2)**2)
	# A_list
	def a(j13, j14):
		t_coord = np.where((coord-j14) % geom['n_phi'] == 0)[0]
		return np.sum(temp[t_coord]*np.exp(-2.0*pi*1j*coord*j13/geom['n_phi']))
	A_list = np.array([a(x[0], x[1]) for x in product(range(1-geom['n_phi'], geom['n_phi']), repeat=2)],dtype=np.complex128)
	A_list = A_list.reshape(2*geom['n_phi']-1, 2*geom['n_phi']-1)
	return A_list.conj() * inter['inten']

def twobodycoulomb_sqrt(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")
	coord = np.arange(1 - up_bound, up_bound)
	d = inter['d']

	q2inter = _twobodyq2(geom, inter['ani_inter'])
	q2mass = _twobodyq2(geom, inter['ani_mass'])
	#q2inter[up_bound-1, up_bound-1] = pi
	temp = np.exp(-q2mass/2)
	temp = temp*(2*pi/np.sqrt(q2inter + inter['d']**2))/(2*geom['area'])
	#temp[up_bound-1, up_bound-1] = 0
	aa = np.zeros(inter['level']+1)
	aa[inter['level']] = 1.0
	_interaction = np.polynomial.laguerre.Laguerre(aa)
	temp *= (_interaction(q2mass/2)**2)
	# A_list
	def a(j13, j14):
		t_coord = np.where((coord-j14) % geom['n_phi'] == 0)[0]
		return np.sum(temp[t_coord]*np.exp(-2.0*pi*1j*coord*j13/geom['n_phi']))
	A_list = np.array([a(x[0], x[1]) for x in product(range(1-geom['n_phi'], geom['n_phi']), repeat=2)],dtype=np.complex128)
	A_list = A_list.reshape(2*geom['n_phi']-1, 2*geom['n_phi']-1)
	return A_list.conj() * inter['inten']

def twobodycoulomb_exp(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")
	coord = np.arange(1 - up_bound, up_bound)
	d = inter['d']

	q2inter = _twobodyq2(geom, inter['ani_inter'])
	q2mass = _twobodyq2(geom, inter['ani_mass'])
	q2inter[up_bound-1, up_bound-1] = pi
	temp = np.exp(-q2mass/2)
	temp = np.exp(-inter['d']*np.sqrt(q2inter))*temp*(2*pi/np.sqrt(q2inter))/(2*geom['area'])
	temp[up_bound-1, up_bound-1] = 0
	aa = np.zeros(inter['level']+1)
	aa[inter['level']] = 1.0
	_interaction = np.polynomial.laguerre.Laguerre(aa)
	temp *= (_interaction(q2mass/2)**2)
	# A_list
	def a(j13, j14):
		t_coord = np.where((coord-j14) % geom['n_phi'] == 0)[0]
		return np.sum(temp[t_coord]*np.exp(-2.0*pi*1j*coord*j13/geom['n_phi']))
	A_list = np.array([a(x[0], x[1]) for x in product(range(1-geom['n_phi'], geom['n_phi']), repeat=2)],dtype=np.complex128)
	A_list = A_list.reshape(2*geom['n_phi']-1, 2*geom['n_phi']-1)
	return A_list.conj() * inter['inten']

def twobodypseudo2body(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")
	coord = np.arange(1 - up_bound, up_bound)

	if inter['order'] < 0:
		coef = np.ones(geom['pseudo2'],dtype=np.float64)
		_interaction = np.polynomial.laguerre.Laguerre(coef*inter['inten'])
	else:
		coef = np.zeros(inter['order']+1,dtype=np.float64)
		coef[inter['order']] = inter['inten']
		_interaction = np.polynomial.laguerre.Laguerre(coef)

	q2 = _twobodyq2(geom, inter['ani_intrin'])
	temp = np.exp(-q2/2)
	temp = temp*(_interaction(q2))/(2*geom['area'])
	
	aa = np.zeros(inter['level']+1)
	aa[inter['level']] = 1.0
	_interaction = np.polynomial.laguerre.Laguerre(aa)
	temp *= (_interaction(q2/2)**2)
	# A_list
	def a(j13, j14):
		t_coord = np.where((coord-j14) % geom['n_phi'] == 0)[0]
		return np.sum(temp[t_coord]*np.exp(-2*pi*1j*coord*j13/geom['n_phi']))
	A_list = np.array([a(x[0], x[1]) for x in product(range(1-geom['n_phi'], geom['n_phi']), repeat=2)],dtype=np.complex128)
	A_list = A_list.reshape(2*geom['n_phi']-1, 2*geom['n_phi']-1)
	return A_list.conj()

def twobodygpseudo2body(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")
	coord = np.arange(1 - up_bound, up_bound)

	_interaction = genlaguerre(inter['orderm'], inter['ordern'])

	q2 = _twobodyq2(geom, inter['ani_intrin'])
	vq = _twobodyvq(geom)
	temp = np.exp(-q2/2)
	temp = temp*(_interaction(q2))*2.0*np.real(vq**inter['ordern'])/(2*geom['area'])
	NN = 2**(inter['ordern']-1)
	for i in range(inter['orderm']+1, inter['orderm']+inter['ordern']+1):
		NN /= i
	NN = np.sqrt( NN/np.pi )
	temp *= NN
	if inter['ordern']==0:
		temp /= np.sqrt(2)

	# A_list
	def a(j13, j14):
		t_coord = np.where((coord-j14) % geom['n_phi'] == 0)[0]
		return np.sum(temp[t_coord]*np.exp(-2*pi*1j*coord*j13/geom['n_phi']))
	A_list = np.array([a(x[0], x[1]) for x in product(range(1-geom['n_phi'], geom['n_phi']), repeat=2)],dtype=np.complex128)
	A_list = A_list.reshape(2*geom['n_phi']-1, 2*geom['n_phi']-1)
	return A_list.conj()

def threebodyfermion(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")

	bound = 2*geom['n_phi']
	coord = np.arange(-bound, bound+1)
	qx = coord/geom['lx']
	qy = coord/geom['ly'] - coord.reshape(2*bound+1, 1)*math.cos(geom['theta'])/geom['lx']
	qx *= 2*pi
	qy *= 2*pi/math.sin(geom['theta'])
	d = bound - 1
	A_list = np.zeros((d,d,d,d),dtype=np.complex128)
	A_list_t = np.zeros((d,d,d,d),dtype=np.complex128)

	Q, phi = inter['ani_intrin'];
	g11= np.cosh(Q)+np.cos(phi)*np.sinh(Q)
	g12= np.sin(phi)*np.sinh(Q)
	g22= np.cosh(Q)-np.cos(phi)*np.sinh(Q) #P.R.B 155140 (2018) Eq.(2)
	fermion.A_Three_Body_0(A_list_t, A_list, qx, qy, bound, 
							geom['n_phi'], geom['num_thread'], g11, g12, g22)
	gc.collect()
	A_list /= (-geom['area']*geom['area'])
	return A_list.conj()

def threebodyboson(inter, geom):
	if not isinstance(inter, dict):
		raise ValueError("inter is not a dictionary instance!")
	if not isinstance(geom, dict):
		raise ValueError("geom is not a dictionary instance!")

	bound = 2*geom['n_phi']
	coord = np.arange(-bound, bound+1)
	qx = coord/geom['lx']
	qy = coord/geom['ly'] - coord.reshape(2*bound+1, 1)*math.cos(geom['theta'])/geom['lx']
	qx *= 2*pi
	qy *= 2*pi/math.sin(geom['theta'])
	d = bound - 1
	A_list = np.zeros((d,d,d,d),dtype=np.complex128)
	A_list_t = np.zeros((d,d,d,d),dtype=np.complex128)

	Q, phi = inter['ani_intrin'];
	g11= np.cosh(Q)+np.cos(phi)*np.sinh(Q)
	g12= np.sin(phi)*np.sinh(Q)
	g22= np.cosh(Q)-np.cos(phi)*np.sinh(Q) #P.R.B 155140 (2018) Eq.(2)
	boson.A_Three_Body_0(A_list_t, A_list, qx, qy, bound, 
							geom['n_phi'], geom['num_thread'], g11, g12, g22)
	gc.collect()
	A_list /= (geom['area']*geom['area'])
	return A_list.conj()