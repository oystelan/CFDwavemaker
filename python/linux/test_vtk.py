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
#path = os.path.join(temp, "libCFDwavemaker_openmp.so")
print(temp)
#lib = CDLL("./comflow_wavemaker.so")

mydll = cdll.LoadLibrary("../../builds/linux64/libCFDwavemaker_vtk_openmp.so")     

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





time = np.arange(10, 30., 0.05)


ww = []
u = []
#u2 = []
#u3 = []
v = []
w = []

x = -795.
y = 0.
z = -250.

for t in time:
    #print(waveelev(mydll,t,x,y))

    ww.append(waveelev(mydll,t,x,y))
    u.append(velocityX(mydll,t,x,y,z))
    #u2.append(velocityX(mydll,t,x,y,z-20.))
    #u3.append(velocityX(mydll,t,x,y,z-40.))
    v.append(velocityY(mydll,t,x,y,z))
    w.append(velocityZ(mydll,t,x,y,z))
    




plt.clf()
plt.plot(time,u, label='u')
plt.plot(time,v, label='v')
plt.plot(time,w, label='w')
#plt.plot(time, ww, label="eta")
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig("./result_u_w_eta.png")


clean_up(mydll)


