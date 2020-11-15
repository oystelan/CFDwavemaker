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
path = os.path.join(temp, "libCFDwavemaker.so")
print(temp)
#lib = CDLL("./comflow_wavemaker.so")

mydll = cdll.LoadLibrary("./libCFDwavemaker.so")     

#exit()

def velocityX(mydll,t,x,y,z):
    aa = mydll.VelocityX
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))
    
def velocityY(mydll,t,x,y,z):
    aa = mydll.VelocityY
    aa.restype = c_double
    aa.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def velocityZ(mydll,t,x,y,z):
    aa = mydll.VelocityZ
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def waveelev(mydll,t,x,y):
    aa = mydll.SurfaceElevation
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))


def volfrac(mydll,x,y,z,t,delta):
    aa = mydll.VolumeFraction
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(delta))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0))

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa()

print(init_dll(mydll))


time = np.arange(0,100,0.5)

z = []
u = []
v = []
w = []
for t in time:
    ww = waveelev(mydll,t,0,0)
    print(t, ww)
    z.append(ww)
    u.append(velocityX(mydll,t,0,0,0))
    v.append(velocityY(mydll,t,0,0,0))
    w.append(velocityZ(mydll,t,0,0,0))




