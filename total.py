#-*-coding:utf-8-*-
from util import *
close("all")
figure(figsize=(8,8))
ax=gca()
p=0.1;s=0.80
ax.set_position([p,p,s,s])

########################################
#MAP LIMITS
########################################
lat0=4
lon0=-75
dlat=30
dlon=30

########################################
#CREATE MAP
########################################
m=Map(projection='stere',
      resolution='l',
      lat_0=lat0,lon_0=lon0,
      llcrnrlon=lon0-dlon,llcrnrlat=lat0-dlat,
      urcrnrlon=lon0+dlon,urcrnrlat=lat0+dlat)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawmapboundary(fill_color='aqua')
m.drawparallels(arange(-80.0,81.0,5.0),labels=[1,1,0,0],fontsize=8)
m.drawmeridians(arange(-180.0,181.0,5.0),labels=[0,0,1,1],fontsize=8)

########################################
#SHOW LIMITS AT A GIVEN LONGITUDE
########################################
alt=0.0;

filename="data/limits-alt_%.5e-lon0_%.4lf-dlon_%.4lf.dat"%(alt,lon0,dlon)
deltal=1
if not fileexists(filename):
    data=[]
    for lon in arange(lon0-dlon,lon0+dlon+deltal,deltal):
        print "Calculating limits for lon = %.4lf..."%lon
        out=System("./contact-latitudes.out %f %f 2> tmp/latitudes.log"%(lon,alt)).strip()
        lats,latn=[float(x) for x in out.split(',')]
        xn,yn=m(lon,latn)
        xs,ys=m(lon,lats)
        data+=[[lon,lats,latn,xn,yn,xs,ys]]
    data=array(data)
    savetxt(filename,data)
else:
    print "File already exists..."
    data=loadtxt(filename)

bn=interp(data[:,3],data[:,4],kind='slinear')
bs=interp(data[:,5],data[:,6],kind='slinear')
xmin=max(data[:,3].min(),data[:,5].min())
xmax=min(data[:,3].max(),data[:,5].max())
xs=linspace(xmin,xmax,100)
fill_between(xs,bs(xs),bn(xs),color='k',alpha=0.4,zorder=+10)

title(r"Ocultacion de Marte por la Luna, Julio 5 de 2014, 21:00 a 22:00 UTC-5",
      position=(0.5,1.05),fontsize=12)
text(0.5,-0.08,r"SAA - UdeA (Medellin - Colombia)",transform=ax.transAxes,
     horizontalalignment='center',fontsize=10)
text(0.5,-0.11,r"Tiempos de contacto en: http://bit.ly/aristarco-0705-contactos",transform=ax.transAxes,
     horizontalalignment='center',fontsize=10)

savefig("images/totalidad.png")
