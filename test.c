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
  double params[100];
  params[0]=lon;
  params[1]=lat;
  params[2]=alt;
  params[3]=0.0;

  double t;

  SpiceDouble dmars,RAmars,DECmars,
    RAmarsJ2000,DECmarsJ2000,ltmars;
  SpiceDouble dmoon,RAmoon,DECmoon,RApole,DECpole,q,
    RAmoonJ2000,DECmoonJ2000,ltmoon;
  SpiceDouble dRA;
  double angdist;
  SpiceDouble PA,PAp;
  SpiceDouble extra[100];

  //////////////////////////////////////////
  //CONTACT TIME
  //////////////////////////////////////////
  //GET TMIN
  int status;double tmin,dcenter;SpiceChar utc[100];
  occultationCord(lon,lat,alt,0,&status,&tmin,&dcenter);
  et2utc_c(tmin+dut*HOUR,"C",2,100,utc);
  printf("tmin = %s\n",utc);

  //FIRST CONTACT
  double t1;SpiceChar utc_t1[100];
  //*
  params[3]=-1.0;
  t1=contactTime(tmin,TEND,params);//CONTACT 3
  //*/
  /*
  params[3]=0.0;
  t1=contactTime(tmin,TEND,params);//CONTACT 3.5
  //*/
  /*
  params[3]=+1.0;
  t1=contactTime(TINI,tmin,params);//CONTACT 1
  //*/
  //*
  params[3]=+1.0;
  t1=contactTime(tmin,TEND,params);//CONTACT 4
  //*/
  et2utc_c(t1+dut*HOUR,"C",2,100,utc_t1);
  fprintf(stderr,"Contacto: %s\n",utc_t1);

  printf("c = %.8e, t1 = %.17e\n",CPARAM,t1);

  //t1=4.57889016099993587e+08;
  t=t1;
  bodyEphemeris("MOON","MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		 &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);

  printf("Moon: RA = %s, DEC = %s, lt = %e\n",dec2sex(RAmoonJ2000),dec2sex(DECmoonJ2000),ltmoon/60);
  printf("Moon: RA = %s, DEC = %s, lt = %e\n",dec2sex(RAmoon),dec2sex(DECmoon),ltmoon/60);
  exit(0);

  //////////////////////////////////////////
  //POSITION ANGLE
  //////////////////////////////////////////
   //////////////////////////////////////////
  //EPHEMERIS INCLUDING POLE POSITION
  //////////////////////////////////////////
  /*
  SpiceDouble extra[10],bp;
  int ie=0;
  extra[ie++]=RMOON;extra[ie++]=FMOON;
  bodyEphemeris2("MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		 &RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon,extra);
  RApole=extra[ie++];DECpole=extra[ie++];q=extra[ie++];
  printf("Moon: Center = (%s,%s), Pole = (%s,%s)\n",
	 dec2sex(RAmoon),dec2sex(DECmoon),dec2sex(RApole),dec2sex(DECpole));
  angdist=R2D(greatCircleDistance(D2R(RAmoon*15),D2R(RApole*15),D2R(DECmoon),D2R(DECpole)));
  printf("Dist: %e (%e)\n",angdist,params[14]);
  PAp=positionAngle(RApole*15,RAmoon*15,DECpole,DECmoon);
  printf("PA pole: %e\n",PAp);
  printf("Pole angle: %e\n",q);
  
  //CALCULATE b'
  bp=BMOON*sqrt(sin(D2R(q))*sin(D2R(q))+(AMOON/BMOON)*(AMOON/BMOON)*cos(D2R(q))*cos(D2R(q)));
  printf("a = %e, b = %e: b' = %e\n",AMOON,BMOON,bp);

  bodyEphemeris("MARS",t,CPARAM,lon,lat,alt,&dmars,&ltmars,
		&RAmarsJ2000,&DECmarsJ2000,&RAmars,&DECmars);
  bodyEphemeris("MOON",t,CPARAM,lon,lat,alt,&dmoon,&ltmoon,
		&RAmoonJ2000,&DECmoonJ2000,&RAmoon,&DECmoon);

  angdist=R2D(greatCircleDistance(D2R(RAmoon*15),D2R(RAmars*15),D2R(DECmoon),D2R(DECmars)));
  dRA=(RAmars-RAmoon)*15;
  
  PA=R2D(asin(sin(D2R(dRA))*cos(D2R(DECmars))/sin(D2R(angdist))));

  printf("\n");
  printf("Difference: %e\n",dRA);
  printf("Mars Declination: %e\n",DECmars);
  printf("Angular distance: %e\n",angdist);
  printf("\n");

  PA=positionAngle(RAmars*15,RAmoon*15,DECmars,DECmoon);

  printf("Mars, Moon: RA = (%s,%s), DEC = (%s,%s)\n",
	 dec2sex(RAmars),dec2sex(RAmoon),
	 dec2sex(DECmars),dec2sex(DECmoon));
  printf("Position angle: %e\n",PA);

  //CALCULATE RMOON
  SpiceDouble Rmoon,dPA;
  dPA=PA-PAp;
  printf("PA - PAp = %e\n",dPA);
  Rmoon=sqrt((AMOON*sin(D2R(dPA)))*(AMOON*sin(D2R(dPA)))+(BMOON*cos(D2R(dPA)))*(BMOON*cos(D2R(dPA))));
  printf("Rmoon = %s, RMOON = %s\n",dec2sex(R2D(Rmoon/dmoon)),dec2sex(R2D(RMOON/dmoon)));
  */
  //////////////////////////////////////////
  //RETURN
  //////////////////////////////////////////
  return 0;
}
