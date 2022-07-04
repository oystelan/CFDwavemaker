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

mydll = cdll.LoadLibrary("../../builds/linux64/libCFDwavemaker_swd_openmp.so")     

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


clean_up(mydll)

exit()

time = np.arange(0,20,0.1)


ww = []
ww2 = []
u = []
v = []
w = []

for t in time:
    print(waveelev(mydll,t,-200,0))

    ww.append(waveelev(mydll,t,-200,0))
    ww2.append(waveelev(mydll,t,-200.,3.))
    u.append(velocityX(mydll,t,-200.,0.,-5.))
    v.append(velocityY(mydll,t,-200.,0.,-5.))
    w.append(velocityZ(mydll,t,-200.,0.,-5.))
    

plt.plot(time, ww, label="x=-200, y=0")
plt.plot(time, ww2, label="x=-150, y=5.")
plt.legend()
plt.grid(True)
plt.savefig("./result_eta.png")

plt.clf()
plt.plot(time,u, label='u')
plt.plot(time,v, label='v')
plt.plot(time,w, label='w')
plt.legend()
plt.grid(True)
plt.savefig("./result_u.png")





