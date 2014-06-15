"""
All projections:
http://matplotlib.org/basemap/users/mapsetup.html
"""
from mpl_toolkits.basemap import Basemap as Map,shiftgrid as Grid
from matplotlib.pyplot import *
from numpy import *
from sys import exit

def drawMap(proj='robin',
            proj_opts=dict(resolution='c',lon_0=0),
            pars=arange(-60.,115,15.),
            pars_opts=dict(labels=[1,1,0,0],labelstyle="+/-",fontsize=8),
            mers=arange(-360.,360.,30.),
            mers_opts=dict(labels=[0,0,1,1],labelstyle="+/-",fontsize=8),
            coasts=False,
            fill=False,
            coasts_opts=dict(linewidth=0.5),
            fill_opts=dict(color='g',alpha=0.3,lake_color='aqua')
            ):

    m=Map(projection=proj,**proj_opts)
    m.drawmapboundary()
    parallels=m.drawparallels(pars,**pars_opts)
    meridians=m.drawmeridians(mers,**mers_opts)
    if coasts:
        m.drawcoastlines(**coasts_opts)
    if fill:
        m.fillcontinents(**fill_opts)
    return m

def plotMap(map,alpha,delta,**args):
    x,y=map(alpha,delta)
    plot(x,y,**args)

def textMap(map,alpha,delta,txt,**args):
    x,y=map(alpha,delta)
    text(x,y,txt,**args)

def readStars():
    import csv
    f=open("stars.csv")
    data=csv.reader(f)
    i=0
    
    stars=[]
    props=[]
    for row in data:
        i+=1
        if i==1:continue

        try:
            BpV=float(row[13])
        except:
            BpV=-111

        #TEXTUAL
        stars+=[
            dict(
                id=row[0],
                hip=row[1],
                hd=row[2],
                hr=row[3],
                gl=row[4],
                denom=row[5],
                name=row[6],
                sp=row[12],
                RA=float(row[7]),
                dec=float(row[8]),
                dist=float(row[9]),
                mag=float(row[10]),
                Mag=float(row[11]),
                BpV=BpV
                )
            ]

        #NUMERIC
        props+=[[
            float(row[7]),#0:RA
            float(row[8]),#1:dec
            float(row[9]),#2:dist
            float(row[10]),#3:mag
            float(row[11]),#4:Mag
            BpV,#5:B-V
            ]]
        i+=1

    return array(stars),array(props)
