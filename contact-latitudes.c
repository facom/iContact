#include <util.c>

#define NARGS 2

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
  double alt=atof(argv[2]);
  fprintf(stderr,"Testing longitude: %.5lf, altitude: %.5e\n",lon,alt);

  //////////////////////////////////////////
  //GET LIMITING LATITUDES
  //////////////////////////////////////////
  double latn,lats,dlat=0.1;
  
  double param[2];
  param[0]=lon;
  param[1]=alt;
  latn=bisectEquation(occultationDistance,param,0.0,90.0,1E-3);
  fprintf(stderr,"Northern Latitude:%.4lf\n",latn);
  lats=bisectEquation(occultationDistance,param,-90.0,0.0,1E-3);
  fprintf(stderr,"Southern Latitude:%.4lf\n",lats);
  fprintf(stderr,"Limits in latitude: (%.4lf,%.4lf)\n",lats,latn);

  //////////////////////////////////////////
  //REPORT RESULTS
  //////////////////////////////////////////
  fprintf(stdout,"%.5lf,%.5lf\n",lats,latn);
  
  return 0;
}
