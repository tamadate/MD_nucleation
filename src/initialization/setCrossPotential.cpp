#include "../variables.hpp"

void
Variables::setCrossPotentials(void){
  int Natypes=atypes.size();
  pair_coeff.resize(Natypes);
  for (int i=0;i<Natypes;i++){
    pair_coeff[i].resize(Natypes);
    for (int j=0;j<Natypes;j++){
      pair_coeff[i][j].resize(2);
    }
  }
  for (int i=0;i<Natypes;i++){
    for (int j=0;j<Natypes;j++){
      double epu=sqrt(atypes[i].coeff1*atypes[j].coeff1);
      double sigma=(atypes[i].coeff2+atypes[j].coeff2)*0.5;
      //sigma=sqrt(atypes_v[i].coeff2*atypes_v[j].coeff2);
      pair_coeff[i][j][0]=48 * epu*pow(sigma,12.0);
      pair_coeff[i][j][1]=24 * epu*pow(sigma,6.0);
    }
  }

}
