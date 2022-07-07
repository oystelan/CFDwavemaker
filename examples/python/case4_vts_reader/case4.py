# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
import os

# give location of dll
#mydll = cdll.comflow
temp = os.path.abspath(__file__)
temp = os.path.realpath(temp)
temp = os.path.dirname(temp)

# we load the shared library 
mydll = cdll.LoadLibrary("../../../builds/linux64/libCFDwavemaker_all_openmp.so")     

# and define the functions which CFDwavemaker.h tells us are in the library

def wave_VeloX(mydll,t,x,y,z):
    aa = mydll.wave_VeloX
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_VeloY(mydll,t,x,y,z):
    aa = mydll.wave_VeloY
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_VeloZ(mydll,t,x,y,z):
    aa = mydll.wave_VeloZ
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_SurfElev(mydll,t,x,y):
    aa = mydll.wave_SurfElev
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(t))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0))

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa()


# All done. lets test:

# we start by initializing (reads data from waveinput.dat) 
init_dll(mydll)

x = 4.8
y = 0.23
z = -10.0


time = np.arange(0,8,0.01)

elev = []
u = []
v = []
w = []

for t in time:
    elev.append(wave_SurfElev(mydll,t,x,y))
    u.append(wave_VeloX(mydll,t,x,y,z))
    v.append(wave_VeloY(mydll,t,x,y,z))
    w.append(wave_VeloZ(mydll,t,x,y,z))
    

# now, lets view the resulting surface elevation
plt.plot(time, elev, label="surface elevation")
plt.legend()
plt.xlabel('Time [sec]')
plt.ylabel('Surface elevation [m]')
plt.grid(True)
plt.savefig("./result_eta.png")

# and also the wave kinematics
plt.clf()
plt.plot(time,u, label='u')
plt.plot(time,v, label='v')
plt.plot(time,w, label='w')
plt.xlabel('Time [sec]')
plt.ylabel('particle velocity [m/s]')
plt.legend()
plt.grid(True)
plt.savefig("./result_uvw.png")

plt.clf()
ug = np.gradient(u, time)
vg = np.gradient(v, time)
wg = np.gradient(w, time)
plt.plot(time,ug, label='ug')
plt.plot(time,vg, label='vg')
plt.plot(time,wg, label='wg')
plt.xlabel('Time [sec]')
plt.ylabel('particle velocity gradient [m/s^2]')
plt.legend()
plt.grid(True)
plt.savefig("./result_uvw_gradient.png")


# all done. remember to clean up after us.
clean_up(mydll)