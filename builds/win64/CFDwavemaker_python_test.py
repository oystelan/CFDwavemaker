# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# give location of dll
mydll = cdll.CFDwavemaker
def tic():
    import time
    #Homemade version of matlab tic and toc functions    
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time    
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")

def velocityx(mydll,t,x,y,z):
    aa = mydll.VelocityX
    aa.restype = c_double
    #aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))    


    
def velocityy(mydll,t,x,y,z):
    aa = mydll.VelocityY
    aa.restype = c_double
    #aa.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)]
    aa.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def velocityz(mydll,t,x,y,z):
    aa = mydll.VelocityZ
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def waveelev(mydll,t,x,y):
    aa = mydll.SurfaceElevation
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_double,c_double,c_double]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))

def volfrac(mydll,x,y,z,t,delta):
    aa = mydll.VolumeFraction
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double,c_double]
    #print aa(byref(c_int(0)),byref(c_int(0)),byref(c_int(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)),byref(c_double(0)))    
    return aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(delta))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0));

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa();

plt.style.use('seaborn-pastel')


fig = plt.figure()
ax = plt.axes(xlim=(-200, 200), ylim=(-6, 6))
line, = ax.plot([], [], lw=3)

print(init_dll(mydll))
time = np.arange(0, 100, 0.5)
ypos = np.arange(-200, 200, 10)

tic()

ww = np.zeros(len(ypos))

def init():
    line.set_data([], [])
    return line,


def animate(i):
    for i_y, y in enumerate(ypos):
        ww[i_y] = waveelev(mydll, time[i], 0, y)

    line.set_data(ypos, ww)
    return line,


toc()



anim = FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)

anim.save('test.mp4', writer='ffmpeg')

##result1= mydll.add(10,1)
##result2= mydll.sub(10,1)
##print "Addition value:-"+result1
##print "Substraction:-"+result2

#print z


clean_up(mydll)

#exit()