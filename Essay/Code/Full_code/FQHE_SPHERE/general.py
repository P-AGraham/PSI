# coding: utf-8
import numpy as np
from scipy import sparse
from numpy import pi
from ctypes import *
from pathlib import Path
import time, math, sys, cmath, gc, os, copy

from .tools.check import ispositiveint

abpath = "/home/hu/d/ubuntu/TEMPFILE/"
#abpath = "/lab_mem/huliangdong/FQHE_ED_DATA/"

def get_abpath():
	return abpath

def set_abpath(newpath):
	abpath = newpath
	return

class FQHE():
	def __init__(self, n, n_o, num_thread=8):
		if ispositiveint(n):
			self.__n = n
		else:
			raise ValueError("The particle number n must be a positive integer!")

		if ispositiveint(n_o):# n_phi=2s0+1
			self.__n_o = n_o
		else:
			raise ValueError("The orbital number n_phi must be a positive integer!")

		if ispositiveint(num_thread):
			self.__num_thread = num_thread
		else:
			raise ValueError("The number of threads must be a positive integer!")

	@property
	def n(self):
		return self.__n

	def reset_n(self, n):
		self.__n = n

	@property
	def n_o(self):
		return self.__n_o

	@property
	def num_thread(self):
		return self.__num_thread

def mkdir(path):
	path = path.strip()
	path = path.rstrip("\\")
	isExists = os.path.exists(path)
	if not isExists:
		os.makedirs(path)
		return True



from .core.Hilbertspace.QHFMFermion.QHFMFermion_wrapper import cg_py, wigner_3j_py

def cg(la, lb, l, m1, m2, m):
	return cg_py(la, lb, l, m1, m2, m)

def wigner_3j(la, lb, l, m1, m2, m):
	return wigner_3j_py(la, lb, l, m1, m2, m)