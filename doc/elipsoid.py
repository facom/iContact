from matplotlib.pyplot import *
from numpy import *
from sys import exit,argv

def Mrot(ty,tz):
    M=[[cos(ty)*cos(tz),-cos(ty)*sin(tz),-sin(ty)],
       [sin(tz),cos(tz),0],
       [sin(ty)*cos(tz),-sin(ty)*sin(tz),cos(ty)]]
    return array(M)

def iMrot(ty,tz):
    M=[[cos(ty)*cos(tz),(-1)*sin(tz),sin(ty)*cos(tz)],
       [(-1)*-cos(ty)*sin(tz),cos(tz),(-1)*-sin(ty)*sin(tz)],
       [-sin(ty),0,cos(ty)]]
    return array(M)

a=1.0
b=0.5

N=50
thetas=linspace(0,pi,N)
phis=linspace(0,2*pi,N)

ty=float(argv[1])*pi/180
#ty=30*pi/180
tz=0*pi/180
M=iMrot(ty,tz)

data=[]
table=[]
data2=[]
proy=[]
for j in xrange(len(phis)):
    rhomax=0
    xmax=0;ymax=0
    for i in xrange(len(thetas)):
        x=a*sin(thetas[i])*cos(phis[j])
        y=a*sin(thetas[i])*sin(phis[j])
        z=b*cos(thetas[i])
        rho=sqrt(x**2+y**2+z**2)
        table+=[[thetas[i],rho]]
        data+=[[x,y,z]]
        vec=dot(M,array([x,y,z]))
        data2+=[vec]
        rho=sqrt(vec[0]**2+vec[1]**2)
        if rho>rhomax:
            rhomax=rho
            xmax=vec[0];ymax=vec[1];zmax=vec[2]
            tmax=thetas[i]
    if j==0 or True:
        proy+=[[xmax,ymax,zmax,tmax]]

q=arctan(a/b*1/tan(ty))
print q*180/pi
proy=array(proy)
bpp=b*sqrt(sin(ty)**2+(a/b)**2*cos(ty)**2)
print bpp
print max(abs(proy[:,0]))
exit(0)
savetxt("projection.dat",proy)
savetxt("elipsoid.dat",data)
savetxt("elipsoid2.dat",data2)
savetxt("distance.dat",table)

mer=[]
j=0
for i in xrange(len(thetas)):
    x=a*sin(thetas[i])*cos(phis[j])
    y=a*sin(thetas[i])*sin(phis[j])
    z=b*cos(thetas[i])
    rho=sqrt(x**2+y**2+z**2)
    mer+=[dot(M,array([x,y,z]))]
savetxt("meridian.dat",mer)

#phis[j]=3*pi/2+0.1
q=arctan(b/a*(1/tan(ty))*(1/cos(phis[j])))
#if phis[j]>pi/2:q=pi+q
#print "Theta = %f"%(q*180/pi)
#exit(0)

fl=open("point.dat","w")
x=a*sin(q)*cos(phis[j])
y=a*sin(q)*sin(phis[j])
z=b*cos(q)
vec=dot(M,array([x,y,z]))
fl.write("%f %f %f\n"%(vec[0],vec[1],vec[2]))
fl.close()

inters=[]
N=200
phis=linspace(0,2*pi,N)
for j in xrange(len(phis)):
    q=arctan(b/a*(1/tan(ty))*(1/cos(phis[j])))
    if phis[j]>pi/2 and phis[j]<3*pi/2:q=pi+q
    x=a*sin(q)*cos(phis[j])
    y=a*sin(q)*sin(phis[j])
    z=b*cos(q)
    inters+=[dot(M,array([x,y,z]))]
inters=array(inters)
#print max(abs(inters[:,0]))
savetxt("intersection.dat",inters)

bp=b/sqrt(sin(ty)**2+(b/a)**2*cos(ty)**2)
alpha=linspace(0,2*pi,100)
elips=[]
trot=tz-pi/2
M2=array([[cos(trot),-sin(trot)],[sin(trot),cos(trot)]])
for i in xrange(len(alpha)):
    x=a*cos(alpha[i])
    y=bp*sin(alpha[i])
    elips+=[dot(M2,[x,y])]
savetxt("elipse.dat",elips)

for ty in linspace(0,80.0,10):
    ty=ty*pi/180.0
    inters=[]
    N=200
    phis=linspace(0,2*pi,N)
    for j in xrange(len(phis)):
        q=arctan(b/a*(1/tan(ty))*(1/cos(phis[j])))
        if phis[j]>pi/2 and phis[j]<3*pi/2:q=pi+q
        x=a*sin(q)*cos(phis[j])
        y=a*sin(q)*sin(phis[j])
        z=b*cos(q)
        inters+=[dot(M,array([x,y,z]))]
    inters=array(inters)
    brot=max(abs(inters[:,0]))
    #print ty*180/pi,brot

plane=[]
fl=open("plane.dat","w")
for x in linspace(-1,1,10):
    for y in linspace(-1,1,10):
        fl.write("%f %f 0\n"%(x,y))
    fl.write("\n")
fl.close()        

zvec=[0,0,0.7]
zvec=dot(M,array(zvec))

fl=open("z.dat","w")
fl.write("0 0 0 %f %f %f\n"%(zvec[0],zvec[1],zvec[2]))
fl.close()

