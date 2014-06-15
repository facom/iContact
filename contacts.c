#include <util.c>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

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

int main(void)
{
  //////////////////////////////////////////
  //CONSTANTS
  //////////////////////////////////////////
  initSpice();

  //////////////////////////////////////////
  //READ LOCATION
  //////////////////////////////////////////
  FILE *fl;
  SpiceDouble lon_d,lon_m,lon_s,lon;
  SpiceDouble lat_d,lat_m,lat_s,lat;
  SpiceDouble alt;
  fl=fopen("location.conf","r");
  fscanf(fl,"%lf %lf %lf",&lon_d,&lon_m,&lon_s);
  lon=SGN(lon_d)*(fabs(lon_d)+lon_m/60.0+lon_s/3600.0);
  fscanf(fl,"%lf %lf %lf",&lat_d,&lat_m,&lat_s);
  lat=SGN(lat_d)*(fabs(lat_d)+lat_m/60.0+lat_s/3600.0);
  fscanf(fl,"%lf",&alt);
  fclose(fl);

  //////////////////////////////////////////
  //TIMES
  //////////////////////////////////////////
  SpiceDouble t,tini,tend;
  str2et_c("07/06/2014 01:30:00.000",&tini);
  tend=tini+2*HOUR;

  //////////////////////////////////////////
  //VERIFY OCCULTATION
  //////////////////////////////////////////
  FILE *fp=fopen("cord.dat","w");
  SpiceDouble dmars,RAmars,DECmars,RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RAmoonJ2000,DECmoonJ2000,ltmoon;
  double dcenter,dcentermin=1E100,tmin,sizemin;
  for(t=tini;t<=tend;t+=1*MINUTE){
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
    fprintf(fp,"%e %e %e %e %e %e %e %e %e\n",tend,dmars,RMARS,RAmars*15,DECmars,dmoon,RMOON,RAmoon*15,DECmoon);
  }
  printf("Angular size of the moon: %e\n",sizemin);
  printf("Minimum moon Center Distance: %e\n",dcentermin);
  printf("Time: %e\n",tmin);
  
  if(dcentermin>sizemin){
    printf("In this location there are no occultation.\n");
    exit(0);
  }
  
  //////////////////////////////////////////
  //CONTACT TIME
  //////////////////////////////////////////
  SpiceChar utc[100];
  double t1,tm1,t2,t3,tm3,t4,duracion;
  double params[4];
  params[0]=lon;
  params[1]=lat;
  params[2]=alt;
  params[3]=0.0;

  FILE *fc=fopen("icontacts.dat","w");
  double ic;
  for(ic=1E-10;ic<=2;ic+=0.1){
    CPARAM=CSPEED/ic;

    params[3]=+1.0;
    t1=contactTime(tini,tmin,params);
    et2utc_c(t1,"C",2,100,utc);
    printf("Contacto 1: %s\n",utc);
    
    params[3]=0.0;
    tm1=contactTime(tini,tmin,params);
    et2utc_c(tm1,"C",2,100,utc);
    printf("Mitad del contacto 1: %s\n",utc);
    
    params[3]=-1.0;
    t2=contactTime(tini,tmin,params);
    et2utc_c(t2,"C",2,100,utc);
    printf("Contacto 2: %s\n",utc);
    
    params[3]=-1.0;
    t3=contactTime(tmin,tend,params);
    et2utc_c(t3,"C",2,100,utc);
    printf("Contacto 3: %s\n",utc);
    
    params[3]=0.0;
    tm3=contactTime(tmin,tend,params);
    et2utc_c(tm3,"C",2,100,utc);
    printf("Mitad del contacto 3: %s\n",utc);
    
    params[3]=+1.0;
    t4=contactTime(tmin,tend,params);
    et2utc_c(t4,"C",2,100,utc);
    printf("Contacto 4: %s\n",utc);

    duracion=tm3-tm1;
    printf("Duracion del Transito: %s\n",dec2sex(duracion/3600.0));

    fprintf(fc,"%e %e %e %e %e %e %e\n",
	    ic,
	    t1-tini,tm1-tini,t2-tini,
	    t3-tini,tm3-tini,t4-tini);

  }
  fclose(fc);

  return 0;
}
