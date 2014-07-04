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
  //RETURN
  //////////////////////////////////////////
  return 0;
}
