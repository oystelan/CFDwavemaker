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
from time import sleep
from scipy.io import savemat, loadmat



# give location of dll
#mydll = cdll.comflow
temp = os.path.abspath(__file__)
temp = os.path.realpath(temp)
temp = os.path.dirname(temp)

dlclose_func = cdll.LoadLibrary('').dlclose
dlclose_func.argtypes = [c_void_p]

def clamp(aa, a1, a2):
    return min(max(aa, a1), a2)
# and define the functions which CFDwavemaker.h tells us are in the library
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

time = np.arange(0,20,0.05)
x = 0.0
y = 0.0
dz = 0.1
z = np.arange(-20., 12, dz)

# # Step 1 - Lest compute the exact second order solution
# # we load the shared library 
# mydll1 = dlopen("../../builds/linux64/libCFDwavemaker_swd_openmp.so")


# shutil.copyfile("waveinput_exact.dat", 'waveinput.dat')
# init_dll(mydll1)

# #mydll1 = reload_lib(mydll1)
# #init_dll(mydll1)


# wave_elev1 = []
# umom1 = []
# vmom1 = []
# wmom1 = []
# for t in time:
#     print('time: ', t)
#     #print(wave_SurfElev(mydll,t,-200,0))
#     wave_elev1.append(wave_SurfElev(mydll1,t,x,y))
#     umom_temp = 0.
#     vmom_temp = 0.
#     wmom_temp = 0.
#     for zz in z:
#         u_temp = wave_VeloX(mydll1,t,x,y,zz)
#         v_temp = wave_VeloY(mydll1,t,x,y,zz)
#         w_temp = wave_VeloZ(mydll1,t,x,y,zz)
#         ratio = clamp((wave_elev1[-1]-(zz-dz/2.)) /((zz+dz/2.)-(zz-dz/2.)),0., 1.)
#         umom_temp += u_temp*abs(u_temp)*dz*ratio
#         vmom_temp += v_temp*abs(v_temp)*dz*ratio
#         wmom_temp += w_temp*abs(w_temp)*dz*ratio
#     umom1.append(umom_temp)
#     vmom1.append(vmom_temp)
#     wmom1.append(wmom_temp)

# clean_up(mydll1)
# dlclose(mydll1)

# data = {}
# data['time'] = time
# data['dz'] = dz
# data['x'] = x
# data['y'] = y
# data['z'] = z
# data['wave_elev1'] = wave_elev1
# data['umom1'] = umom1
# data['vmom1'] = vmom1
# data['wmom1'] = wmom1
# savemat('exact_solution.mat',data)

#exit()
# instead of calculating (save time), load precompute solution:
data = loadmat('exact_solution.mat')
#time = data['time']
#dz = data['dz']
#x = data['x']
#y = data['y']
#z = data['z']
wave_elev1 = np.squeeze(data['wave_elev1'])
umom1 = np.squeeze(data['umom1'])
vmom1 = np.squeeze(data['vmom1'])
wmom1 = np.squeeze(data['wmom1'])


#sleep(1.)

# Step 2 - now lest do the LSgrid with linear interpolation
mydll2 = dlopen("../../builds/linux64/libCFDwavemaker_swd_openmp.so")     
shutil.copyfile("waveinput_linear.dat", 'waveinput.dat')
init_dll(mydll2)

wave_elev2 = []
umom2 = []
vmom2 = []
wmom2 = []
for t in time:
    print('time: ', t)
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev2.append(wave_SurfElev(mydll2,t,x,y))
    umom_temp = 0.
    vmom_temp = 0.
    wmom_temp = 0.
    for zz in z:
        u_temp = wave_VeloX(mydll2,t,x,y,zz)
        v_temp = wave_VeloY(mydll2,t,x,y,zz)
        w_temp = wave_VeloZ(mydll2,t,x,y,zz)
        ratio = clamp((wave_elev2[-1]-(zz-dz/2.)) /((zz+dz/2.)-(zz-dz/2.)),0., 1.)
        umom_temp += u_temp*abs(u_temp)*dz*ratio
        vmom_temp += v_temp*abs(v_temp)*dz*ratio
        wmom_temp += w_temp*abs(w_temp)*dz*ratio
    umom2.append(umom_temp)
    vmom2.append(vmom_temp)
    wmom2.append(wmom_temp)
clean_up(mydll2)
dlclose(mydll2)
#sleep(1)

print("test")
# Step 3 - finally lets do the LSgrid with spline interpolation
mydll3 = cdll.LoadLibrary("../../builds/linux64/libCFDwavemaker_swd_openmp.so")
#mydll3 = dlopen() 
print("test")
shutil.copyfile("waveinput_spline.dat", 'waveinput.dat')
print("test")
init_dll(mydll3)

wave_elev3 = []
umom3 = []
vmom3 = []
wmom3 = []
for t in time:
    print('time: ', t)
    #print(wave_SurfElev(mydll,t,-200,0))
    wave_elev3.append(wave_SurfElev(mydll3,t,x,y))
    umom_temp = 0.
    vmom_temp = 0.
    wmom_temp = 0.
    for zz in z:
        u_temp = wave_VeloX(mydll3,t,x,y,zz)
        v_temp = wave_VeloY(mydll3,t,x,y,zz)
        w_temp = wave_VeloZ(mydll3,t,x,y,zz)
        ratio = clamp((wave_elev3[-1]-(zz-dz/2.)) /((zz+dz/2.)-(zz-dz/2.)),0., 1.)
        umom_temp += u_temp*abs(u_temp)*dz*ratio
        vmom_temp += v_temp*abs(v_temp)*dz*ratio
        wmom_temp += w_temp*abs(w_temp)*dz*ratio
    umom3.append(umom_temp)
    vmom3.append(vmom_temp)
    wmom3.append(wmom_temp)
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
plt.plot(time,umom1, label='u exact')
plt.plot(time,vmom1, label='v exact')
plt.plot(time,wmom1, label='w exact')
plt.plot(time,umom2,'--', label='u linear')
plt.plot(time,vmom2,'--', label='v linear')
plt.plot(time,wmom2,'--', label='w linear')
plt.plot(time,umom3,':', label='u spline')
plt.plot(time,vmom3,':', label='v spline')
plt.plot(time,wmom3,':', label='w spline')
plt.legend()
plt.grid(True)
plt.savefig("./result_umom.png")
plt.show()
# all done. remember to clean up after us.






