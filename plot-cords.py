from util import *
from sys import argv

myid=argv[1]

fig=figure(figsize=(8,8))
ax=fig.add_axes([0.0,0.0,1.0,1.0],axisbg='k')
#ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')

rmin=1E100
rmax=-1E100
for cord in argv[2:]:
    data=loadtxt(cord)

    #TIMES
    times=data[:,0]
    
    #MARS
    dms=data[:,1]
    Rms=data[:,2]
    RAms=data[:,3]
    DECms=data[:,4]
    
    #MOON
    dMs=data[:,5]
    RMs=data[:,6]
    RAMs=data[:,7]
    DECMs=data[:,8]

    #RELATIVE POSITIONS
    dRAs=RAms-RAMs
    dDECs=DECms-DECMs

    #ANGULAR SIZES
    sizeM=arctan(RMs[0]/dMs[0])*180/pi
    sizem=arctan(Rms[0]/dms[0])*180/pi

    rmin=min(rmin,-1.5*sizeM,min(dRAs),min(dDECs))
    rmax=max(rmax,1.5*sizeM,max(dRAs),max(dDECs))

    plot(dDECs,dRAs,'y-',linewidth=1)

    xM,yM=genCircle((0,0),sizeM)
    plot(yM,xM,'w-',linewidth=1)

"""
rmin=-1.5*sizeM
rmax=+1.5*sizeM
"""

#DIRECTION OF MOTION
hw=0.02;hl=0.03

ipos=0
arrow(dDECs[ipos],dRAs[ipos],
      -(dDECs[ipos]-dDECs[ipos+1]),
      -(dRAs[ipos]-dRAs[ipos+1]),
      head_width=hw,head_length=hl,fc='y',ec='y')
ipos=int(len(dDECs)/2)
arrow(dDECs[ipos],dRAs[ipos],
      -(dDECs[ipos]-dDECs[ipos+1]),
      -(dRAs[ipos]-dRAs[ipos+1]),
      head_width=hw,head_length=hl,fc='y',ec='y')
ipos=-2
arrow(dDECs[ipos],dRAs[ipos],
      -(dDECs[ipos]-dDECs[ipos+1]),
      -(dRAs[ipos]-dRAs[ipos+1]),
      head_width=hw,head_length=hl,fc='y',ec='y')

text(0.95*dDECs[-1],0.95*dRAs[-1],"Marte",color='y')

infoname=cord.replace("cord-","info-")
data=loadtxt(infoname)
title("Cuerda para el sitio long. %+.4f, lat. %+.4f, Alt. %.2f"%(data[0],data[1],data[2]),
      position=(0.5,0.95),color='w')

#ORIENTATION
arrow(0.0,0.0,0.0,-0.1,head_width=hw,head_length=hl,fc='w',ec='w')
text(-0.02,-0.12,"Oeste",color='w')
arrow(0.0,0.0,0.1,0.0,head_width=hw,head_length=hl,fc='w',ec='w')
text(0.12,0.05,"Norte",color='w')

xlim((rmax,rmin))
ylim((rmax,rmin))
grid(color='w')
savefig("tmp/cords-%s.png"%myid)
