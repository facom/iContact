from matplotlib import use
use('Agg')
from mpl_toolkits.basemap import Basemap as Map,shiftgrid as Grid
from matplotlib.pyplot import *
from matplotlib.patches import *
from numpy import *
import commands
from os import system
from os.path import lexists as fileexists
from scipy.interpolate import interp1d as interp
from sys import exit

cosd=lambda x:cos(x*pi/180)
sind=lambda x:sin(x*pi/180)
arccosd=lambda x:arccos(x)*180/pi
arcsind=lambda x:arcsin(x)*180/pi

def System(cmd,out=True):
    if not out:
        system(cmd)
        output=""
    else:
        output=commands.getoutput(cmd)
    return output

def genCircle(center,radius,npoints=1000):
    """
    npoints should be always even
    """
    #CENTER
    ao=center[0]
    fo=center[1]

    #RANGE OF PHI
    fmin=fo-radius
    fmax=fo+radius

    fs=linspace(fmin,fmax,npoints/2+1)

    #VALUES OF ALPHA
    fvec=zeros(npoints)
    avec=zeros(npoints)
    i=0
    for f in fs:
        arg=(cosd(radius)-sind(f)*sind(fo))/(cosd(f)*cosd(fo))
        if arg>1:arg=1.0
        if arg<-1:arg=-1.0
        da=arccosd(arg)
        fvec[i]=f
        avec[i]=ao+da
        fvec[-(i+1)]=f
        avec[-(i+1)]=ao-da
        i+=1

    return avec,fvec

def distanceSphere(center,x,y):
    r=arccosd(sind(y)*sind(center[1])+
              cosd(y)*cosd(center[1])*cosd(x-center[0]))
    return r

