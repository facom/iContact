/*
Kernels:
http://naif.jpl.nasa.gov/pub/naif/
*/

//////////////////////////////////////////
//HEADERS
//////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <SpiceUsr.h>

//////////////////////////////////////////
//CONSTANTS
//////////////////////////////////////////
//WGS84
#define C_REARTH 6378.137
#define C_FEARTH (1/298.257223563)
#define D2R(x) (x*M_PI/180)
#define R2D(x) (x*180/M_PI)
#define POWI(x,n) gsl_pow_int(x,n)

#define SGN(x) (x<0?-1:+1)
#define HOUR 3600.0
#define MINUTE 60.0
#define CSPEED (clight_c())

//////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////
double CPARAM;
double REARTH;
double FEARTH;
double RMOON;
double RMARS;

//////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////
int initSpice(void)
{
  SpiceInt n;
  SpiceDouble radii[3];

  //KERNELS
  furnsh_c("kernels.txt");

  //INITIALIZE GLOBAL VARIABLES
  CPARAM=clight_c();

  //EARTH RADII
  bodvrd_c("EARTH","RADII",3,&n,radii);
  REARTH=radii[0];
  FEARTH=(radii[0]-radii[2])/radii[0];

  //MARS RADII
  bodvrd_c("MARS","RADII",3,&n,radii);
  RMARS=radii[0];

  //MOON RADII
  bodvrd_c("MOON","RADII",3,&n,radii);
  RMOON=radii[0];

  return 0;
}

char *dec2sex(double dec)
{
  double d,m,s;
  int id,im,sgn;
  char *str=calloc(sizeof(char),100); 
  d=fabs(dec);
  sgn=dec/d;
  id=floor(d);
  m=(d-id)*60;
  im=floor(m);
  s=(m-im)*60;
  sprintf(str,"%+d:%02d:%.3f",sgn*id,im,s);
  return str;
}

char *vec2str(double vec[])
{
  char frm[]="%.8e";
  char format[100];
  char *str=calloc(sizeof(char),100); 
  sprintf(format,"%s %s %s",frm,frm,frm);
  sprintf(str,format,vec[0],vec[1],vec[2]);
  return str;
}

int bodyEphemeris(ConstSpiceChar *body,
		  SpiceDouble t,
		  SpiceDouble cspeed,
		  SpiceDouble lon,SpiceDouble lat,SpiceDouble alt,
		  SpiceDouble *range,
		  SpiceDouble *ltime,
		  SpiceDouble *raJ2000,
		  SpiceDouble *decJ2000,
		  SpiceDouble *ra,
		  SpiceDouble *dec
		  )
{
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble bodyJ2000[6],bodySSBJ2000[6],ltbody;
  SpiceDouble bodyTOPOJ2000[3],bodyTOPOEpoch[3];
  SpiceDouble Dbody,RAbody,DECbody,RAbodyJ2000,DECbodyJ2000;
  SpiceDouble observerITRF93[3],observerJ2000[3],observerSSBJ2000[3];
  SpiceDouble M_J2000_Epoch[3][3];
  SpiceDouble M_ITRF93_J2000[3][3];

  SpiceDouble d,lt,ltmp,ltold,lttol=1E-2;
  int i,ncn=10;

  //ROTATION MATRIX AT THE TIME OF EPHEMERIS
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  pxform_c("IAU_EARTH","J2000",t,M_ITRF93_J2000);

  //OBSERVER POSITION J2000
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);

  //LIGHT TIME CORRECTED POSITION
  i=0;
  lt=0.0;ltold=1.0;
  while((fabs(lt-ltold)/lt)>=lttol && i<ncn){
    ltold=lt;
    spkezr_c("EARTH",t,"J2000","NONE","SOLAR SYSTEM BARYCENTER",earthSSBJ2000,&ltmp);
    spkezr_c(body,t-lt,"J2000","NONE","SOLAR SYSTEM BARYCENTER",bodySSBJ2000,&ltmp);
    vsub_c(bodySSBJ2000,earthSSBJ2000,bodyJ2000);
    /*
      //IT SHOULD BE RESPECT TO OBSERVER.  NASA'S HORIZONS IS WRONG
      vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
      vsub_c(bodySSBJ2000,observerSSBJ2000,bodyJ2000);
    //*/
    d=vnorm_c(bodyJ2000);
    lt=d/cspeed;
    i++;
  }

  //TOPOCENTRIC POSITION
  vsub_c(bodyJ2000,observerJ2000,bodyTOPOJ2000);
  //RA & DEC J2000
  recrad_c(bodyTOPOJ2000,&d,&RAbodyJ2000,&DECbodyJ2000);
  //PRECESS POSITION
  mxv_c(M_J2000_Epoch,bodyTOPOJ2000,bodyTOPOEpoch);
  //RA & DEC PRECESSED
  recrad_c(bodyTOPOEpoch,&d,&RAbody,&DECbody);

  *range=d;
  *ltime=lt;
  *ra=RAbody*180/M_PI/15;
  *dec=DECbody*180/M_PI;
  *raJ2000=RAbodyJ2000*180/M_PI/15;
  *decJ2000=DECbodyJ2000*180/M_PI;
  return 0;
}

double greatCircleDistance(double lam1,double lam2,
			   double phi1,double phi2)
{
  double d;

  d=acos(sin(phi1)*sin(phi2)+
	 cos(phi1)*cos(phi2)*cos(lam2-lam1));
  
  return d;
}

double angularRadius(double R,double d)
{
  double aR;

  aR=atan(R/d);

  return aR;
}
