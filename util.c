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

#define EARTH_ID "EARTH"
#define MOON_ID "MOON"
#define MARS_ID "MARS BARYCENTER"
#define MARS_CENTER_ID "MARS"

#define LON 0
#define LAT 1
#define ALT 2
#define TRUE 1
#define FALSE 0

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
double RMOON,AMOON,BMOON,FMOON;
double RMARS,AMARS,BMARS,FMARS;
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
  bodvrd_c(EARTH_ID,"RADII",3,&n,radii);
  REARTH=radii[0];
  FEARTH=(radii[0]-radii[2])/radii[0];

  //MARS RADII
  bodvrd_c(MARS_CENTER_ID,"RADII",3,&n,radii);
  AMARS=radii[0];
  BMARS=radii[2];
  FMARS=(AMARS-BMARS)/AMARS;
  RMARS=(AMARS+BMARS)/2;

  //MOON RADII
  bodvrd_c(MOON_ID,"RADII",3,&n,radii);
  BMOON=1737.10;
  AMOON=1738.14;
  FMOON=(AMOON-BMOON)/AMOON;
  RMOON=(AMOON+BMOON)/2;

  /*
  printf("f: earth = %e, moon = %e, mars = %e\n",
	 FEARTH,FMOON,FMARS);
  */

  //DATE AND TIME OF OCCULTATION
  str2et_c("07/06/2014 01:30:00.000",&TINI);
  TEND=TINI+3*HOUR;

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
		   ConstSpiceChar *bodyname,
		   SpiceDouble t,
		   SpiceDouble cspeed,
		   SpiceDouble lon,SpiceDouble lat,SpiceDouble alt,
		   SpiceDouble *range,
		   SpiceDouble *ltime,
		   SpiceDouble *raJ2000,
		   SpiceDouble *decJ2000,
		   SpiceDouble *ra,
		   SpiceDouble *dec,
		   SpiceDouble *extra
		  )
{
  SpiceDouble earthSSBJ2000[6];
  SpiceDouble bodyJ2000[6],bodySSBJ2000[6],ltbody;
  SpiceDouble bodyTOPOJ2000[3],bodyTOPOEpoch[3];
  SpiceDouble Dbody,RAbody,DECbody,RAbodyJ2000,DECbodyJ2000;
  SpiceDouble observerITRF93[3],observerJ2000[3],observerSSBJ2000[3];
  SpiceDouble M_J2000_Epoch[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  SpiceDouble M_ITRF93_J2000[3][3];

  SpiceDouble d,lt,ltmp,ltold,lttol=1E-2;
  int i,ie=0,ncn=10;

  //ROTATION MATRIX AT THE TIME OF EPHEMERIS
  pxform_c("J2000","EARTHTRUEEPOCH",t,M_J2000_Epoch);
  //pxform_c("IAU_EARTH","J2000",t,M_ITRF93_J2000);
  pxform_c("ITRF93","J2000",t,M_ITRF93_J2000);

  //OBSERVER POSITION J2000
  georec_c(D2R(lon),D2R(lat),alt/1000.0,REARTH,FEARTH,observerITRF93);
  mxv_c(M_ITRF93_J2000,observerITRF93,observerJ2000);

  //LIGHT TIME CORRECTED POSITION
  i=0;
  lt=0.0;ltold=1.0;
  spkezr_c(EARTH_ID,t,"J2000","NONE","SOLAR SYSTEM BARYCENTER",earthSSBJ2000,&ltmp);
  vadd_c(earthSSBJ2000,observerJ2000,observerSSBJ2000);
  while((fabs(lt-ltold)/lt)>=lttol && i<ncn){
    ltold=lt;
    spkezr_c(body,t-lt,"J2000","NONE","SOLAR SYSTEM BARYCENTER",bodySSBJ2000,&ltmp);
    vsub_c(bodySSBJ2000,observerSSBJ2000,bodyTOPOJ2000);
    d=vnorm_c(bodyTOPOJ2000);
    lt=d/cspeed;
    i++;
  }
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

  //POLE POSITION
  SpiceDouble Rbody,fbody,poleBODY[3],poleJ2000[3],M_MOON_J2000[3][3],bodyPoleJ2000[3];
  SpiceDouble bodyPoleTOPOJ2000[3],bodyPoleTOPOEpoch[3];
  SpiceDouble RApole,DECpole;
  SpiceChar Bref[100];
  Rbody=extra[ie++];fbody=extra[ie++];
  georec_c(0.0,M_PI/2,0.0,Rbody,fbody,poleBODY);
  sprintf(Bref,"IAU_%s",bodyname);
  pxform_c(Bref,"J2000",t-lt,M_MOON_J2000);
  mxv_c(M_MOON_J2000,poleBODY,poleJ2000);
  vadd_c(bodyJ2000,poleJ2000,bodyPoleJ2000);
  //printf("POLE: %s\nPOLE J2000: %s\nBODY: %s\nPOLE: %s\n",vec2str(poleBODY),vec2str(poleJ2000),vec2str(bodyJ2000),vec2str(bodyPoleJ2000));

  //VISUAL-POLE ANGLE
  SpiceDouble ub[3],nb,up[3],np,cosq,q;
  unorm_c(bodyJ2000,ub,&nb);
  unorm_c(poleJ2000,up,&np);
  cosq=vdot_c(ub,up);
  q=R2D(acos(cosq));
  if(q>90) q=180-q;

  //POLE COORDINATES
  vsub_c(bodyPoleJ2000,observerJ2000,bodyPoleTOPOJ2000);
  mxv_c(M_J2000_Epoch,bodyPoleTOPOJ2000,bodyPoleTOPOEpoch);
  recrad_c(bodyPoleTOPOEpoch,&d,&RApole,&DECpole);

  extra[ie++]=RApole*180/M_PI/15;
  extra[ie++]=DECpole*180/M_PI;
  extra[ie++]=q;
  return 0;
}

double greatCircleDistance(double lam1,double lam2,
			   double phi1,double phi2)
{
  double d;

  //COSINE FORMULA
  /*
  d=acos(sin(phi1)*sin(phi2)+
	 cos(phi1)*cos(phi2)*cos(lam2-lam1));
  //*/
  
  //HARVESINE FORMULA
  double sf,sl;
  sf=sin((phi2-phi1)/2);
  sl=sin((lam2-lam1)/2);
  d=2*asin(sqrt(sf*sf+cos(phi1)*cos(phi2)*sl*sl));

  return d;
}

double positionAngle(double lam1,double lam2,
		     double phi1,double phi2)
{
  double d,dlam,PA;
  d=R2D(greatCircleDistance(D2R(lam1),D2R(lam2),
			    D2R(phi1),D2R(phi2)));
  dlam=fabs(lam1-lam2);
  PA=R2D(asin(sin(D2R(dlam))*cos(D2R(phi1))/sin(D2R(d))));

  if(phi1<phi2) PA=180-PA;
  if(lam2>lam1) PA=360-PA;

  return PA;
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
  double extra[100];
  int ie=0;
  /*k: Contact parameter
    k = +0: Centers
    k = +1: Outer contact
    k = -1: Inner contact
   */
  double k=ps[3];
  
  extra[0]=RMARS;extra[1]=FMARS;
  bodyEphemeris(MARS_ID,MARS_CENTER_ID,t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		&RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars,extra);
  extra[0]=RMOON;extra[1]=FMOON;
  bodyEphemeris(MOON_ID,MOON_ID,t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		&RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);

  angdist=R2D(greatCircleDistance(D2R(RAmoon*15),D2R(RAmars*15),D2R(DECmoon),D2R(DECmars)));
  //printf("Angular distance: %.17e\n",angdist);exit(0);

  //RADIUS ASSUMING NON-SPHERICITY
  SpiceDouble RApole,DECpole,q,PA,PAp,dPA,bp;
  RApole=extra[ie++];DECpole=extra[ie++];q=extra[ie++];
  PAp=positionAngle(RApole*15,RAmoon*15,DECpole,DECmoon);
  bp=BMOON*sqrt(sin(D2R(q))*sin(D2R(q))+(AMOON/BMOON)*(AMOON/BMOON)*cos(D2R(q))*cos(D2R(q)));
  PA=positionAngle(RAmars*15,RAmoon*15,DECmars,DECmoon);
  dPA=PA-PAp;
  aRM=R2D(angularRadius(sqrt((AMOON*sin(D2R(dPA)))*(AMOON*sin(D2R(dPA)))+(BMOON*cos(D2R(dPA)))*(BMOON*cos(D2R(dPA)))),dmoon));

  //aRM=R2D(angularRadius(RMOON,dmoon));
  aRm=R2D(angularRadius(RMARS,dmars));
  cfunc=angdist-aRM-k*aRm;

  //RETURNING INFO
  int ip=4;
  ps[ip++]=RAmars;//4
  ps[ip++]=DECmars;//5
  ps[ip++]=RAmarsJ2000;//6
  ps[ip++]=DECmarsJ2000;//7
  ps[ip++]=dmars;//8
  ps[ip++]=ltmars;//9
  ps[ip++]=RAmoon;//10
  ps[ip++]=DECmoon;//11
  ps[ip++]=RAmoonJ2000;//12
  ps[ip++]=DECmoonJ2000;//13
  ps[ip++]=dmoon;//14
  ps[ip++]=ltmoon;//15
  ps[ip++]=aRM;//16
  ps[ip++]=aRm;//17
  ps[ip++]=cfunc;//18

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
    status=gsl_root_test_residual(cfunc,1E-6);
    if(status==GSL_SUCCESS) break;
  }while(status==GSL_CONTINUE && niter<maxiter);

  /*
  printf("Moon coordinates (J2000):\n\tRA = %s, DEC = %s, d = %.10e, lt = %.5e, Ang.Size = %.3f\n",
	 dec2sex(params[12]),dec2sex(params[13]),params[14],params[15]/60,2*params[16]*3600);
  printf("Moon coordinates (Epoch):\n\tRA = %s, DEC = %s\n",
	 dec2sex(params[10]),dec2sex(params[11]));
  printf("Mars coordinates (J2000):\n\tRA = %s, DEC = %s, d = %.10e, lt = %.5e, Ang.Size = %.3f\n",
	 dec2sex(params[6]),dec2sex(params[7]),params[8],params[9]/60,2*params[17]*3600);
  printf("Mars coordinates (Epoch):\n\tRA = %s, DEC = %s\n",
	 dec2sex(params[4]),dec2sex(params[5]));
  printf("Cfunc = %.17e\n",params[18]);
  //*/
  
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
  double extra[100];
  int ie=0;
  
  for(t=TINI;t<=TEND;t+=1*MINUTE){
    extra[0]=RMARS;extra[1]=FMARS;
    bodyEphemeris(MARS_ID,"MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		  &RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars,extra);
    extra[0]=RMOON;extra[1]=FMOON;
    bodyEphemeris(MOON_ID,"MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		  &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);
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
  double extra[100];
  int ie=0;
  
  if(verbose){
    sprintf(filename,"data/cord-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    sprintf(infoname,"data/info-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    fp=fopen(filename,"w");
    fi=fopen(infoname,"w");
  }
  for(t=TINI;t<=TEND;t+=1*MINUTE){
    extra[0]=RMARS;extra[0]=FMARS;
    bodyEphemeris(MARS_ID,"MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		  &RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars,extra);
    extra[0]=RMOON;extra[0]=FMOON;
    bodyEphemeris(MOON_ID,"MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		  &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);
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

int occultationCord2(double *location,double tini,double tend,double dt,
		     int verbose,
		     int *Status,double *Tmin,double *Dcenter)
{
  char filename[300],infoname[300];
  FILE *fp,*fi;
  SpiceDouble dmars,RAmars,DECmars,RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RAmoonJ2000,DECmoonJ2000,ltmoon;
  double dcenter,dcentermin=1E100,t,tmin,sizemin;
  double lon=location[LON],lat=location[LAT],alt=location[ALT];
  double extra[100];
  int ie=0;
  
  if(verbose){
    sprintf(filename,"data/cord-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    sprintf(infoname,"data/info-%+.4lf_%+.4lf_%+.3lf.dat",lon,lat,alt);
    fp=fopen(filename,"w");
    fi=fopen(infoname,"w");
  }
  for(t=tini;t<=tend;t+=dt){
    extra[0]=RMARS;extra[0]=FMARS;
    bodyEphemeris(MARS_ID,"MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		  &RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars,extra);
    extra[0]=RMOON;extra[0]=FMOON;
    bodyEphemeris(MOON_ID,"MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		  &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);
    dcenter=R2D(greatCircleDistance(D2R(RAmars)*15,D2R(RAmoon)*15,D2R(DECmars),D2R(DECmoon)));
    if(dcenter<=dcentermin){
      tmin=t;
      dcentermin=dcenter;
      sizemin=R2D(atan(RMOON/dmoon));
    }
    if(verbose) 
      fprintf(fp,"%.17e %e %e %e %e %e %e %e %e %e\n",t,
	      dmars,RMARS,RAmars*15,DECmars,
	      dmoon,RMOON,RAmoon*15,DECmoon,dcenter);
  }
  if(verbose){
    fprintf(fi,"%e %e %e\n",lat,lon,alt);
    fclose(fi);
    fclose(fp);
  }

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
