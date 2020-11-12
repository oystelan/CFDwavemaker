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

def velocityx(mydll,t,x,y,z):
    aa = mydll.VelocityX
    aa.restype = c_double
    aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))    

def velocityy(mydll,t,x,y,z):
    aa = mydll.VelocityY
    aa.restype = c_double
    aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t)) 

def velocityz(mydll,t,x,y,z):
    aa = mydll.VelocityZ
    aa.restype = c_double
    aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t)) 
    
def waveelev(mydll,t,x,y):
    aa = mydll.SurfaceElevation
    aa.restype = c_double
    aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))    

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

#print(waveelev(mydll,5.,0,0))


exit()
time = np.arange(0,10,0.05)

z = []
u = []
v = []
w = []
for t in time:
    z.append(waveelev(mydll,t,0,0))
    u.append(velocityx(mydll,t,0,0,0))
    v.append(velocityy(mydll,t,0,0,0))
    w.append(velocityz(mydll,t,0,0,0))

#plt.plot(time,z,time,u,time,v,time,w)
#plt.plot(time,z)
#plt.show()
##result1= mydll.add(10,1)
##result2= mydll.sub(10,1)
##print "Addition value:-"+result1
##print "Substraction:-"+result2


