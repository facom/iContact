#include <util.c>

#define NARGS 4

int main(int argc,char *argv[])
{
  //////////////////////////////////////////
  //INITIALIZE SPICE
  //////////////////////////////////////////
  initSpice();

  //////////////////////////////////////////
  //VARIABLES
  //////////////////////////////////////////
  int status;
  double loc1[100],loc2[100];
  double tini,tend,deltat,t,tmin,dcenter;

  //////////////////////////////////////////
  //LOCATIONS
  //////////////////////////////////////////
  //LOCATION 1
  loc1[LON]=-75.0;
  loc1[LAT]=1;
  loc1[ALT]=0;

  //LOCATION 2
  loc2[LON]=-75.0;
  loc2[LAT]=-18;
  loc2[ALT]=0;

  //////////////////////////////////////////
  //TIMES
  //////////////////////////////////////////
  str2et_c("07/06/2014 01:30:00.000",&tini);
  tend=tini+3*HOUR;
  deltat=10*MINUTE;

  occultationCord2(loc1,tini,tend,deltat,FALSE,&status,&tmin,&dcenter);
  loc1[3]=+1.0;
  t=contactTime(tini,tmin,loc1);

  tini=t;
  tend=tini+1.5*HOUR;
  
  //////////////////////////////////////////
  //CORDS
  //////////////////////////////////////////
  occultationCord2(loc1,tini,tend,deltat,TRUE,&status,&tmin,&dcenter);
  occultationCord2(loc2,tini,tend,deltat,TRUE,&status,&tmin,&dcenter);

  //////////////////////////////////////////
  //READ DATA
  //////////////////////////////////////////
  FILE *fl;
  char file1[300],file2[300];

  sprintf(file1,"data/cord-%+.4lf_%+.4lf_%+.3lf.dat",
	  loc1[LON],loc1[LAT],loc1[ALT]);
  sprintf(file2,"data/cord-%+.4lf_%+.4lf_%+.3lf.dat",
	  loc2[LON],loc2[LAT],loc2[ALT]);

  double tmp,t2;
  double dmars1,RAmars1,DECmars1,dmoon1,RAmoon1,DECmoon1;
  double dmars12,RAmars12,DECmars12,dmoon12,RAmoon12,DECmoon12;
  double dmars2,RAmars2,DECmars2,dmoon2,RAmoon2,DECmoon2;

  fl=fopen(file1,"r");
  fscanf(fl,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&t,
	 &dmars1,&tmp,&RAmars1,&DECmars1,&dmoon1,&tmp,&RAmoon1,&DECmoon1,
	 &tmp);
  fscanf(fl,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&t2,
	 &dmars12,&tmp,&RAmars12,&DECmars12,&dmoon12,&tmp,&RAmoon12,&DECmoon12,
	 &tmp);
  printf("dmoon1 = %.8e\n",dmoon1);
  printf("ltmoon1 = %.8e\n",dmoon1/CPARAM);
  printf("ltmars1 = %.8e\n",dmars1/CPARAM/60);
  printf("RAmars1 = %s\n",dec2sex(RAmars1));
  printf("DECmars1 = %s\n",dec2sex(DECmars1));
  fclose(fl);

  double drmoon,drmars;
  drmars=R2D(greatCircleDistance(D2R(RAmars1),D2R(RAmars12),
				 D2R(DECmars1),D2R(DECmars12)))/(t2-t);
  drmoon=R2D(greatCircleDistance(D2R(RAmoon1),D2R(RAmoon12),
				 D2R(DECmoon1),D2R(DECmoon12)))/(t2-t);

  printf("Observing time = %e\n",t2-t);
  printf("Rate mars = %e\n",drmars*60*3600);
  printf("Rate moon = %e\n",drmoon*3600);

  fl=fopen(file2,"r");
  fscanf(fl,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&t,
	 &dmars2,&tmp,&RAmars2,&DECmars2,&dmoon2,&tmp,&RAmoon2,&DECmoon2,
	 &tmp);
  printf("dmoon2 = %.8e\n",dmoon2);
  printf("RAmars2 = %s\n",dec2sex(RAmars2));
  printf("DECmars2 = %s\n",dec2sex(DECmars2));
  fclose(fl);

  double d=fabs(dmoon2-dmoon1);
  printf("dmoon1 - dmoon2 = %.8e\n",d/dmoon1);

  //////////////////////////////////////////
  //GEOCENTRIC MOON DISTANCE
  //////////////////////////////////////////
  double lt;
  SpiceDouble rmoon[6],D;
  spkezr_c(MOON_ID,t,"J2000","NONE","EARTH",rmoon,&lt);
  D=vnorm_c(rmoon);
  printf("Theoretical distance: %.8e\n",D);
  
  //////////////////////////////////////////
  //PARALLAX
  //////////////////////////////////////////
  double pi;
  pi=R2D(greatCircleDistance(D2R(RAmars1-RAmoon1),D2R(RAmars2-RAmoon2),
			     D2R(DECmars1-DECmoon1),D2R(DECmars2-DECmoon2)));
  printf("Parallax = %e\n",pi);
  //*
  pi=R2D(greatCircleDistance(D2R(RAmoon1),D2R(RAmoon2),
			     D2R(DECmoon1),D2R(DECmoon2)));
  //*/
  //pi=0.0795;
  printf("Parallax = %e\n",pi);

  //UNITARY VECTOR
  double umoon1[3],umoon2[3];
  radrec_c(1.0,D2R(RAmoon1),D2R(DECmoon1),umoon1);
  radrec_c(1.0,D2R(RAmoon2),D2R(DECmoon2),umoon2);

  //////////////////////////////////////////
  //EARTH
  //////////////////////////////////////////
  double Delta;
  Delta=R2D(greatCircleDistance(D2R(loc1[LON]),D2R(loc2[LON]),D2R(loc1[LAT]),D2R(loc2[LAT])));
  printf("Location separation = %e\n",Delta);
  
  SpiceDouble B;
  SpiceDouble obs1[3],obs2[3],dobs[3],dobsJ[3],E2J[3][3];

  pxform_c("ITRF93","J2000",t,E2J);
  georec_c(D2R(loc1[LON]),D2R(loc1[LAT]),loc1[ALT]/1000.0,REARTH,FEARTH,obs1);
  georec_c(D2R(loc2[LON]),D2R(loc2[LAT]),loc2[ALT]/1000.0,REARTH,FEARTH,obs2);
  vsub_c(obs1,obs2,dobs);
  mxv_c(E2J,dobs,dobsJ);

  double uobs[3];
  unorm_c(dobsJ,uobs,&tmp);
  printf("B = %e km\n",tmp);

  B=vnorm_c(dobsJ);
  printf("B = %e km\n",B);

  //////////////////////////////////////////
  //MOON DISTANCE
  //////////////////////////////////////////
  double Do,dDo;
  Do=B/(2*sin(D2R(pi/2)));
  dDo=REARTH/(2*Do)*100;
  printf("Do = %e +/- %.1f%% km\n",Do,dDo);
  printf("Do = %e - %e km\n",Do*(1-dDo/100),Do*(1+dDo/100));

  double Dobs1,Dobs2;
  Dobs1=Do+REARTH;
  Dobs2=sqrt(Do*Do-REARTH*REARTH);
  printf("Dobs = %e - %e km\n",Dobs1,Dobs2);

  double q1,q2,D1,D2;
  q1=acos(vdot_c(uobs,umoon1));
  q2=acos(vdot_c(uobs,umoon2));
  printf("q1 = %e\n",R2D(q1));
  printf("q2 = %e\n",R2D(q2));

  D1=B/sin(D2R(pi))*sin(q2);
  printf("D1 = %e km +/- %.4f %%\n",D1,fabs(D1-dmoon1)/dmoon1*100);
  D2=B/sin(D2R(pi))*sin(q1);
  printf("D2 = %e km +/- %.4f %%\n",D2,fabs(D2-dmoon2)/dmoon2*100);

  return 0;
}
