# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
import os
import shutil



# give location of dll
#mydll = cdll.comflow
temp = os.path.abspath(__file__)
temp = os.path.realpath(temp)
temp = os.path.dirname(temp)

dlclose_func = cdll.LoadLibrary('').dlclose
dlclose_func.argtypes = [c_void_p]

def dlopen(fpath):
    mydll = cdll.LoadLibrary(fpath)
    return mydll

def is_loaded(lib):
    libp = os.path.abspath(lib)
    ret = os.system("lsof -p %d | grep %s > /dev/null" % (os.getpid(), libp))
    return (ret == 0)

def reload_lib(lib):
    handle = lib._handle
    name = lib._name
    del lib
    while is_loaded(name):
        dlclose_func(handle)   
        #libdl = CDLL("libdl.so")
        #libdl.dlclose(handle)
    return cdll.LoadLibrary(name)

# def dlclose(mydll):
#     handle = mydll._handle
#     del mydll
#     dlclose_func(handle)

def dlclose(lib):
    handle = lib._handle
    name = lib._name
    del lib
    while is_loaded(name):
        dlclose_func(handle)   
        #libdl = CDLL("libdl.so")
        #libdl.dlclose(handle)       

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

def wave_Kinematics(mydll,t,x,y,z):
    aa = mydll.wave_Kinematics
    #aa.restype = c_double
    aa.restype = POINTER(c_double*4)
    aa.argtypes = [c_double,c_double,c_double,c_double]
    res = aa(c_double(x),c_double(y),c_double(z),c_double(t)).contents
    #print(res[0], res[1], res[2], res[3])
    return res[0], res[1], res[2], res[3]

def wave_update(mydll, t):
    aa = mydll.wave_force_update
    aa.restype = c_int
    aa.argtype = c_double
    return(aa(c_double(t)))

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

time = np.arange(0,20,0.1)
x = 0.0
y = 0.0
z = -10

# Step 1 - Lest compute the exact second order solution
shutil.copyfile("waveinput_exact.dat", 'waveinput.dat')
mydll1 = dlopen("../../builds/linux64/libCFDwavemaker_swd_openmp.so")  
init_dll(mydll1)

wave_elev1 = []
u1 = []
v1 = []
w1 = []

for t in time:
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev1.append(wave_SurfElev(mydll1,t,x,y))
    u1.append(wave_VeloX(mydll1,t,x,y,z))
    v1.append(wave_VeloY(mydll1,t,x,y,z))
    w1.append(wave_VeloZ(mydll1,t,x,y,z))
clean_up(mydll1)
dlclose(mydll1)

# Step 2 - now lest do the LSgrid with linear interpolation
mydll2 = dlopen("../../builds/linux64/libCFDwavemaker_swd_openmp.so")  
shutil.copyfile("waveinput_linear.dat", 'waveinput.dat')
init_dll(mydll2)

wave_elev2 = []
u2 = []
v2 = []
w2 = []

for t in time:
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev2.append(wave_SurfElev(mydll2,t,x,y))
    u2.append(wave_VeloX(mydll2,t,x,y,z))
    v2.append(wave_VeloY(mydll2,t,x,y,z))
    w2.append(wave_VeloZ(mydll2,t,x,y,z))
clean_up(mydll2)
dlclose(mydll2)

# Step 3 - finally lets do the LSgrid with spline interpolation
mydll3 = dlopen("../../builds/linux64/libCFDwavemaker_swd_openmp.so")  
shutil.copyfile("waveinput_spline.dat", 'waveinput.dat')
init_dll(mydll3)

wave_elev3 = []
u3 = []
v3 = []
w3 = []

for t in time:
    wave_update(mydll3, t)
    data = wave_Kinematics(mydll3,t,x,y,z)
    #print(wave_SurfElev(mydll,t,-200,0))
    # wave_elev3.append(wave_SurfElev(mydll3,t,x,y))
    # u3.append(wave_VeloX(mydll3,t,x,y,z))
    # v3.append(wave_VeloY(mydll3,t,x,y,z))
    # w3.append(wave_VeloZ(mydll3,t,x,y,z))
    wave_elev3.append(data[0])
    u3.append(data[1])
    v3.append(data[2])
    w3.append(data[3])
clean_up(mydll3)
dlclose(mydll3)
    




plt.plot(time, wave_elev1, label="exact")
plt.plot(time, wave_elev2,"--", label="linear")
plt.plot(time, wave_elev3,"--", label="spline")
plt.legend()
plt.grid(True)
plt.savefig("./result_eta.png")
plt.show()
plt.clf()
plt.plot(time,u1, label='u exact')
plt.plot(time,v1, label='v exact')
plt.plot(time,w1, label='w exact')
plt.plot(time,u2,'--', label='u linear')
plt.plot(time,v2,'--', label='v linear')
plt.plot(time,w2,'--', label='w linear')
plt.plot(time,u3,':', label='u spline')
plt.plot(time,v3,':', label='v spline')
plt.plot(time,w3,':', label='w spline')
plt.legend()
plt.grid(True)
plt.savefig("./result_u.png")
plt.show()
# all done. remember to clean up after us.






