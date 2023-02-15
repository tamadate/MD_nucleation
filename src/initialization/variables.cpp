#include "../variables.hpp"

Variables::Variables(void) {
  
  for (int nth=0;nth<Nth;nth++){
    Potentials Us;
    Us.Uion=Us.Ugas=Us.Uvap=Us.Ugi=Us.Ugg=Us.Uvg=Us.Uvi=Us.Uvv=0;
  }
}

void
Variables::ionInitial(double T) {
  random_device seed;
  default_random_engine engine(seed());
  double A,B,C,x,y,z;
  A=seed(),B=seed(),C=seed();
  for(auto &r : positionAtom[0]) {
    x=r[0],y=r[1],z=r[2]; 
    ROTATION(r[0],r[1],r[2],A,B,C,x,y,z);
  }

  int Nion=positionAtom[0].size();
  for(int i=0; i<Nion; i++) {
    double matom=massAtom[0][i]*1e-3/Nw;
    normal_distribution<> dist(0.0, sqrt(kb*T/matom));
    velocityAtom[0][i][0]=dist(engine)*1e-5;
    velocityAtom[0][i][1]=dist(engine)*1e-5;
    velocityAtom[0][i][2]=dist(engine)*1e-5;
  }
}
