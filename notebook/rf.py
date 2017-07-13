
# imports {{{1
import numpy as np
from math import pi
from numpy import fft, random
from functools import reduce, partial

import sys
from heapq import *
#}}}

# box properties {{{1
class Box:
	def __init__(self, dim, N, L):
		self.N = N
		self.L = L
		self.dim = dim
		
		self.shape = (self.N,) * dim
		self.size = reduce(lambda x, y: x*y, self.shape)
	
		self.KN = wave_number(self.shape) * 2*pi / self.N
		self.K = wave_number(self.shape) * 2*pi / self.L
		self.k = np.sqrt((self.K**2).sum(axis=0))

		self.k1 = self.k.copy()
		self.k1[0,0] = 1.0 # the 0.0 here often gives errors
	
		self.k_max = N*pi/L
		self.k_min = 2*pi/L

		self.DX = np.array([[1, 0], [0, 1], [1, -1], [1, 1]])

	def index(self, f, x):
		return f[x[1] % self.shape[1], x[0] % self.shape[0]]

	def print_status(self, f):
		f.write("%f Mpc/h, %u^%u box ============\n" % (self.L, self.N, self.dim))
		f.write("Largest scale: min_k = 2*pi / L = %f h/Mpc\n" % (2*pi / self.L))
		f.write("Smallest scale: max_k = N*pi / L = %f h/Mpc\n" % (self.N*pi / self.L))

#}}}

# random fields {{{1
# we need an array to give coordinates of k-space
def wave_number(s):
	N = s[0]
	i = np.indices(s)
	return np.where(i > N/2, i - N, i)

#def print_array(out, A):
#	for l in A:
#		print >> out, " ".join(map(str, l))
#	print >> out, ''
#	print >> out, ''

# discrete scale-space filter
def scale_filter(K, t, gamma):
	return np.exp(-(2 - gamma) * t + \
		(1 - gamma) * (np.cos(K[0]) + np.cos(K[1])) * \
		t + gamma * np.cos(K[0]) * np.cos(K[1]) * t)

class ScaleFreePS:
	def __init__(self, n):
		self.n = n

	def __call__(self, k):
		power = k**self.n
		power.flat[0] = 0
		return power

def PowerSpectrum(type, args):
#	if type == "BBKS":
#		return BBKS_PS(*args)
#
#	if type == "EH":
#		return EH_PS(*args)

	if type == "scale-free":
		return ScaleFreePS(*args)
#
#	if type == "table":
#		return TablePS(*args)

	print("Unrecognised power spectrum type: %s" % type)
	sys.exit(2)


def random_field(B, n, scale):
	A_white = random.normal(0.0, 1.0, B.shape)
	F = fft.fftn(A_white) 
	
	# multiply the fourier transform with sqrt(P(k))
	P = PowerSpectrum("scale-free", (n ,))
	F *= np.where(B.k > B.k_max, 0, np.sqrt(P(B.k1)))
	F[0,0] = 0

	# scale if needed
	if (scale != None):
		F *= scale_filter(B.KN, scale / B.L * B.N, 4./9)

	# inverse transform
	A = fft.ifftn(F).real
	return A
#}}}

# eigen values {{{1
def Potential(f, B):
	F = fft.fft2(f)
	P = fft.ifft2(- F / B.k1**2).real
	return P

def Hessian(f, B):
	F = fft.fft2(f)
	H11 = fft.ifft2(- F * B.K[0]**2).real
	H22 = fft.ifft2(- F * B.K[1]**2).real
	H12 = fft.ifft2(- F * B.K[1]*B.K[0]).real
	return {11: H11, 12: H12, 22: H22}

def EigenValues(H):
	d = (H[11] + H[22]) / 2.0
	q = np.sqrt((H[11] - H[22])**2 + H[12]**2) / 2.0
	return d + q, d - q
#}}}

