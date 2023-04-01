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

# int bose_shooter(int d, double *params, double xst, double xend, double omega_st, double omega_end, double *params_out)
bose_shooter_cpp = lib.bose_shooter
bose_shooter_cpp.argtypes = [ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'), 
                             ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                             ctl.ndpointer(np.float64, flags='aligned, c_contiguous')]
bose_shooter_cpp.restype = ctypes.c_int

# int get_star_profile(int d, double *params, double xst, double xend, int npoints, double *xarr, double *yarr)
get_star_profile_cpp = lib.get_star_profile
get_star_profile_cpp.argtypes = [ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                                 ctypes.c_double, ctypes.c_double, ctypes.c_int,
                                 ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                                 ctl.ndpointer(np.float64, flags='aligned, c_contiguous')]
get_star_profile_cpp.restype = ctypes.c_int

class Bose_Star_Solver:
    def __init__(self, d = 3, l = 1.0, kappa = 1.0, el = 0.1, mass = 0.0, phi0 = 0.1, xrange = [0.001, 1.4]):
        self.d = d
        self.params = np.array([l, kappa, mass, el, phi0])
        self.xrange = xrange
        
    def shoot_for_omega(self, omega_range = [1, 2]):
        output = np.zeros(5)
        code = bose_shooter_cpp(self.d, self.params, self.xrange[0], self.xrange[1], omega_range[0], omega_range[1], output)
        if(code != 0):
            sys.stdout.flush()
            raise RuntimeError('No omega obtained')
        
        self.omega = output[0]
        self.delta0 = output[1]
        self.a0 = output[2]
        self.star_mass = output[3];
        self.star_charge = output[4];

        print("Omega found by shooting %.4lf" % self.omega)
        if(self.params[1] != 0.0):
            print("Star mass is equal to %.4lf" % self.star_mass)
        if(self.params[3] != 0.0):
            print("Star charge is equal to %.4lf" % self.star_charge)

    def get_star_profile(self, npoints=501):
        self.coords = np.linspace(0.0, np.pi/2.0, npoints)
        self.fields = np.zeros(6*npoints)

        setup_parametes = np.r_[self.params, [self.omega, self.delta0, self.a0]]
        code = get_star_profile_cpp(self.d, setup_parametes, self.xrange[0], self.xrange[1], npoints, self.coords, self.fields)
        if(code != 0):
            raise RuntimeError('No star profile was obtained')

        self.fields = self.fields.reshape(-1,6)

    def plot_star_profile(self):
        fig = plt.figure(figsize=(12,8))

        ax = plt.subplot(221)
        ax.set_title("$\\phi(x)$")
        ax.plot(self.coords, self.fields[:,2])

        ax = plt.subplot(222)
        ax.set_title("$A(x)$")
        ax.plot(self.coords, self.fields[:,0])

        ax = plt.subplot(223)
        ax.set_title("$\\delta(x)$")
        ax.plot(self.coords, self.fields[:,1])

        ax = plt.subplot(224)
        ax.set_title("$a(x)$")
        ax.plot(self.coords, self.fields[:,4])

# int get_difference_with_asymptotic(int d, double *params, double xst, double xend, double *dy_diff)
get_difference_with_asymptotic_cpp = lib.get_difference_with_asymptotic
get_difference_with_asymptotic_cpp.argtypes = [ctypes.c_int, ctl.ndpointer(np.float64, flags='aligned, c_contiguous'),
                                 ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double)]
get_difference_with_asymptotic_cpp.restype = ctypes.c_int

class Bose_Star_Manual:
    def __init__(self, d = 3, l = 1.0, kappa = 1.0, el = 0.1, mass = 0.0, phi0 = 0.1, xrange = [0.001, 1.4]):
        self.d = d
        self.params = np.array([l, kappa, mass, el, phi0])
        self.xrange = xrange

        self.omega = 0.0
        self.delta0 = 0.0
        self.a0 = 0.0

    def shoot_dy(self, omega, delta0=0.0, a0=0.0):
        output = ctypes.c_double(0.0)

        self.omega = omega
        self.delta0 = delta0
        self.a0 = a0

        setup_parametes = np.r_[self.params, [self.omega, self.delta0, self.a0]]

        code = get_difference_with_asymptotic_cpp(self.d, setup_parametes, self.xrange[0], self.xrange[1], output)
        if(code != 0):
            sys.stdout.flush()
            raise RuntimeError('Integration fail')

        return output.value

    def get_star_profile(self, npoints=501):
        self.coords = np.linspace(0.0, np.pi/2.0, npoints)
        self.fields = np.zeros(6*npoints)

        setup_parametes = np.r_[self.params, [self.omega, self.delta0, self.a0]]
        code = get_star_profile_cpp(self.d, setup_parametes, self.xrange[0], self.xrange[1], npoints, self.coords, self.fields)
        if(code != 0):
            raise RuntimeError('No star profile was obtained')

        self.fields = self.fields.reshape(-1,6)

    def plot_star_profile(self):
        fig = plt.figure(figsize=(12,8))

        ax = plt.subplot(221)
        ax.set_title("$\\phi(x)$")
        ax.plot(self.coords, self.fields[:,2])

        ax = plt.subplot(222)
        ax.set_title("$A(x)$")
        ax.plot(self.coords, self.fields[:,0])

        ax = plt.subplot(223)
        ax.set_title("$\\delta(x)$")
        ax.plot(self.coords, self.fields[:,1])

        ax = plt.subplot(224)
        ax.set_title("$a(x)$")
        ax.plot(self.coords, self.fields[:,4])