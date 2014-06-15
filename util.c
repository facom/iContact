/*
Kernels:
http://naif.jpl.nasa.gov/pub/naif/
*/

//////////////////////////////////////////
//HEADERS
//////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <SpiceUsr.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

//////////////////////////////////////////
//MACROS
//////////////////////////////////////////
#define D2R(x) (x*M_PI/180)
#define R2D(x) (x*180/M_PI)
#define POWI(x,n) gsl_pow_int(x,n)
#define SGN(x) (x<0?-1:+1)

//////////////////////////////////////////
//CONSTANTS
//////////////////////////////////////////
//WGS84
#define C_REARTH 6378.137
#define C_FEARTH (1/298.257223563)
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
double TINI,TEND;

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

  //DATE AND TIME OF OCCULTATION
  str2et_c("07/06/2014 01:30:00.000",&TINI);
  TEND=TINI+2*HOUR;

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

double contactFunction(double t,void *params)
{
  SpiceDouble dmars,RAmars,DECmars,RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RAmoonJ2000,DECmoonJ2000,ltmoon;
  double angdist,aRM,aRm;
  double cfunc;
  double* ps=(double*)params;

  double lon=ps[0];
  double lat=ps[1];
  double alt=ps[2];
  /*k: Contact parameter
    k = +0: Centers
    k = +1: Outer contact
    k = -1: Inner contact
   */
  double k=ps[3];
  
  bodyEphemeris("MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		&RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars);
  bodyEphemeris("MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		&RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon);
  angdist=R2D(greatCircleDistance(D2R(RAmoon*15),D2R(RAmars*15),D2R(DECmoon),D2R(DECmars)));
  aRM=R2D(angularRadius(RMOON,dmoon));
  aRm=R2D(angularRadius(RMARS,dmars));

  cfunc=angdist-aRM-k*aRm;
  return cfunc;
}

double contactTime(double tini,double tend,double *params)
{
   //////////////////////////////////////////
  //PREPARE SOLUTION
  //////////////////////////////////////////
  double t;
  int niter,maxiter=100,status;
  double cfunc;
  gsl_root_fsolver *solver;
  solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_function F;
  F.function=&contactFunction;
  F.params=params;

  //////////////////////////////////////////
  //SOLUTION
  //////////////////////////////////////////
  gsl_root_fsolver_set(solver,&F,tini,tend);
  niter=0;
  do{
    niter++;
    status=gsl_root_fsolver_iterate(solver);
    t=gsl_root_fsolver_root(solver);
    cfunc=contactFunction(t,params);
    status=gsl_root_test_residual(cfunc,1E-5);
    if(status==GSL_SUCCESS) break;
  }while(status==GSL_CONTINUE && niter<maxiter);
  
  return t;
}

double occultationDistance(double lat,void *param)
{
  double* ps=(double*)param;
  double lon=ps[0];
  double alt=ps[1];

  SpiceDouble dmars,RAmars,DECmars,RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RAmoonJ2000,DECmoonJ2000,ltmoon;
  double dcenter,dcentermin=1E100,t,tmin,sizemin;
  
  for(t=TINI;t<=TEND;t+=1*MINUTE){
    bodyEphemeris("MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		  &RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars);
    bodyEphemeris("MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		  &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon);
    dcenter=R2D(greatCircleDistance(D2R(RAmars)*15,D2R(RAmoon)*15,D2R(DECmars),D2R(DECmoon)));
    if(dcenter<=dcentermin){
      tmin=t;
      dcentermin=dcenter;
      sizemin=R2D(atan(RMOON/dmoon));
    }
  }
  return (dcentermin-sizemin);
}

int occultationCord(double lon,double lat,double alt,int verbose,
		    int *Status,double *Tmin,double *Dcenter)
{
  char filename[300],infoname[300];
  FILE *fp,*fi;
  SpiceDouble dmars,RAmars,DECmars,RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RAmoonJ2000,DECmoonJ2000,ltmoon;
  double dcenter,dcentermin=1E100,t,tmin,sizemin;
  
  if(verbose){
    sprintf(filename,"data/cord-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    sprintf(infoname,"data/info-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    fp=fopen(filename,"w");
    fi=fopen(infoname,"w");
  }
  for(t=TINI;t<=TEND;t+=1*MINUTE){
    bodyEphemeris("MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		  &RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars);
    bodyEphemeris("MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		  &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon);
    dcenter=R2D(greatCircleDistance(D2R(RAmars)*15,D2R(RAmoon)*15,D2R(DECmars),D2R(DECmoon)));
    if(dcenter<=dcentermin){
      tmin=t;
      dcentermin=dcenter;
      sizemin=R2D(atan(RMOON/dmoon));
    }
    if(verbose) fprintf(fp,"%.17e %e %e %e %e %e %e %e %e %e\n",t,dmars,RMARS,RAmars*15,DECmars,dmoon,RMOON,RAmoon*15,DECmoon,dcenter);
  }
  if(verbose){
    fprintf(fi,"%e %e %e\n",lat,lon,alt);
    fclose(fi);
    fclose(fp);
  }

  /*
  if(verbose){
    printf("Angular size of the moon: %e\n",sizemin);
    printf("Minimum moon Center Distance: %e\n",dcentermin);
    printf("Time: %e\n",tmin);
  }
  //*/
  
  if(dcentermin>sizemin)
    *Status=0;
  else
    *Status=1;

  *Tmin=tmin;
  *Dcenter=dcentermin;

  return 0;
}

double bisectEquation(double (*func)(double,void*),double* params,double pini,double pend,double functol)
{
   //////////////////////////////////////////
  //PREPARE SOLUTION
  //////////////////////////////////////////
  int niter,maxiter=100,status;
  double p,cfunc;
  gsl_root_fsolver *solver;
  solver=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_function F;
  F.function=func;
  F.params=params;

  //////////////////////////////////////////
  //SOLUTION
  //////////////////////////////////////////
  gsl_root_fsolver_set(solver,&F,pini,pend);
  niter=0;
  do{
    niter++;
    status=gsl_root_fsolver_iterate(solver);
    p=gsl_root_fsolver_root(solver);
    cfunc=func(p,params);
    status=gsl_root_test_residual(cfunc,functol);
    if(status==GSL_SUCCESS) break;
  }while(status==GSL_CONTINUE && niter<maxiter);
  
  return p;
}
