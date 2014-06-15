#include <util.c>

#define NARGS 4

int main(int argc,char *argv[])
{
  //////////////////////////////////////////
  //INITIALIZE SPICE
  //////////////////////////////////////////
  initSpice();
  
  //////////////////////////////////////////
  //READ COMMAND LINE ARGUMENTS
  //////////////////////////////////////////
  if(argc<1+NARGS){
    fprintf(stderr,"Number of arguments %d is less than expected (%d).\n",argc-1,NARGS);
    return 1;
  }
  double lon=atof(argv[1]);
  double lat=atof(argv[2]);
  double alt=atof(argv[3]);
  double dut=atof(argv[4]);
  fprintf(stderr,"Contact times for longitude: %.5lf, latitude: %.5lf, altitude: %.5e, UTC%.1lf\n",lon,lat,alt,dut);

  //////////////////////////////////////////
  //GETTING CONTACT TIMES
  //////////////////////////////////////////
  SpiceChar utc1[100],utcm1[100],utc2[100],utc3[100],utcm3[100],utc4[100];
  SpiceChar c_utc1[100],c_utcm1[100],c_utc2[100],c_utc3[100],c_utcm3[100],c_utc4[100];
  SpiceChar dur[100];
  int status;
  double tmin,dcenter;
  double t1,tm1,t2,t3,tm3,t4,duracion;
  double params[4];
  params[0]=lon;
  params[1]=lat;
  params[2]=alt;
  params[3]=0.0;

  //GET TMIN
  occultationCord(lon,lat,alt,0,&status,&tmin,&dcenter);

  //GET CONTACT TIMES
  params[3]=+1.0;
  t1=contactTime(TINI,tmin,params);
  et2utc_c(t1+dut*HOUR,"C",2,100,utc1);
  fprintf(stderr,"Contacto 1: %s\n",utc1);
  
  params[3]=0.0;
  tm1=contactTime(TINI,tmin,params);
  et2utc_c(tm1+dut*HOUR,"C",2,100,utcm1);
  fprintf(stderr,"Mitad del contacto 1: %s\n",utcm1);
  
  params[3]=-1.0;
  t2=contactTime(TINI,tmin,params);
  et2utc_c(t2+dut*HOUR,"C",2,100,utc2);
  fprintf(stderr,"Contacto 2: %s\n",utc2);
  
  params[3]=-1.0;
  t3=contactTime(tmin,TEND,params);
  et2utc_c(t3+dut*HOUR,"C",2,100,utc3);
  fprintf(stderr,"Contacto 3: %s\n",utc3);
  
  params[3]=0.0;
  tm3=contactTime(tmin,TEND,params);
  et2utc_c(tm3+dut*HOUR,"C",2,100,utcm3);
  fprintf(stderr,"Mitad del contacto 3: %s\n",utcm3);
  
  params[3]=+1.0;
  t4=contactTime(tmin,TEND,params);
  et2utc_c(t4+dut*HOUR,"C",2,100,utc4);
  fprintf(stderr,"Contacto 4: %s\n",utc4);
  
  duracion=tm3-tm1;
  strcpy(dur,dec2sex(duracion/3600.0));
  fprintf(stderr,"Duracion del Transito: %s\n",dur);

  //PLOT CORD
  double tbackini=TINI,tbackend=TEND;
  TINI=t1-0.5*HOUR;
  TEND=t4+0.5*HOUR;
  occultationCord(lon,lat,alt,1,&status,&tmin,&dcenter);

  //INFINITE SPEED OF LIGHT
  TINI=tbackini;TEND=tbackend;
  CPARAM=3E10;
  params[3]=+1.0;
  t1=contactTime(TINI,tmin,params);
  et2utc_c(t1+dut*HOUR,"C",2,100,c_utc1);
  fprintf(stderr,"Contacto 1: %s\n",c_utc1);
  
  params[3]=0.0;
  tm1=contactTime(TINI,tmin,params);
  et2utc_c(tm1+dut*HOUR,"C",2,100,c_utcm1);
  fprintf(stderr,"Mitad del contacto 1: %s\n",c_utcm1);
  
  params[3]=-1.0;
  t2=contactTime(TINI,tmin,params);
  et2utc_c(t2+dut*HOUR,"C",2,100,c_utc2);
  fprintf(stderr,"Contacto 2: %s\n",c_utc2);
  
  params[3]=-1.0;
  t3=contactTime(tmin,TEND,params);
  et2utc_c(t3+dut*HOUR,"C",2,100,c_utc3);
  fprintf(stderr,"Contacto 3: %s\n",c_utc3);
  
  params[3]=0.0;
  tm3=contactTime(tmin,TEND,params);
  et2utc_c(tm3+dut*HOUR,"C",2,100,c_utcm3);
  fprintf(stderr,"Mitad del contacto 3: %s\n",c_utcm3);
  
  params[3]=+1.0;
  t4=contactTime(tmin,TEND,params);
  et2utc_c(t4+dut*HOUR,"C",2,100,c_utc4);
  fprintf(stderr,"Contacto 4: %s\n",c_utc4);
  
  //////////////////////////////////////////
  //REPORT RESULTS
  //////////////////////////////////////////
  fprintf(stdout,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
	  utc1,utcm1,utc2,utc3,utcm3,utc4,dur,
	  c_utc1,c_utcm1,c_utc2,c_utc3,c_utcm3,c_utc4);
  
  return 0;
}
