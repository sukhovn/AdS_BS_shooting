import ctypes
import numpy.ctypeslib as ctl
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

if os.path.exists(os.environ["STAR_ROUTINES_PATH"] + '/shooter_lib.so') == False:
	raise ImportError('No shooter_lib.so library found')
else:
	lib=ctl.load_library('shooter_lib.so', os.environ["STAR_ROUTINES_PATH"])

bose_shooter_cpp = lib.bose_shooter
bose_shooter_cpp.argtypes = [ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'), 
                             ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                             ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
bose_shooter_cpp.restype = ctypes.c_int

get_star_profile_cpp = lib.get_star_profile
get_star_profile_cpp.argtypes = [ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'), 
                             ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                             ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                             ctl.ndpointer(np.float64, flags='aligned, c_contiguous')]
get_star_profile_cpp.restype = ctypes.c_int

class Bose_Star_Solver:
    def __init__(self, d = 3, l = 1.0, kappa = 1.0, mass = 0.0, phi0 = 0.1, xrange = [0.001, 1.4]):
        self.d = d
        self.params = np.array([l, kappa, mass, phi0])
        self.xrange = xrange
        
    def shoot_for_omega(self, omega_range = [1, 2]):
        ct_omega = ctypes.c_double(1.0)
        ct_delta0 = ctypes.c_double(0.0)
        
        code = bose_shooter_cpp(self.d, self.params, self.xrange[0], self.xrange[1], omega_range[0], omega_range[1], ct_omega, ct_delta0)
        if(code != 0):
            sys.stdout.flush()
            print("No omega was obtained")
            return
        
        self.omega = ct_omega.value
        self.delta0 = ct_delta0.value
        print("Omega found by shooting %.4lf" % self.omega)
    
    def get_star_profile(self, npoints=501):
        self.coords = np.linspace(0.0, np.pi/2.0, npoints)
        self.fields = np.zeros(4*npoints)
        
        code = get_star_profile_cpp(self.d, self.params, self.xrange[0], self.xrange[1], self.omega, self.delta0, npoints, self.coords, self.fields)
        if(code != 0):
            print("No star profile was obtained")
            return
        
        self.fields = self.fields.reshape(-1,4)
        
    def plot_star_profile(self):
        fig = plt.figure(figsize=(18,4))
        
        ax = plt.subplot(131)
        ax.set_title("$\\phi(x)$")
        ax.plot(self.coords, self.fields[:,2])
        
        ax = plt.subplot(132)
        ax.set_title("$A(x)$")
        ax.plot(self.coords, self.fields[:,0])
                
        ax = plt.subplot(133)
        ax.set_title("$\\delta(x)$")
        ax.plot(self.coords, self.fields[:,1])
        
        plt.show()