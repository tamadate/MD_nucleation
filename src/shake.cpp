#include "md.hpp"
#define tolerance 1e-20
#define rho_zero 0.01

void
MD::updatePositionSHAKE(void){
    vars->times.tpos-=omp_get_wtime();
    int Nmol=vars->mass.size();
    for(int i=0; i<Nmol; i++){
        if(vars->region[i]){
            // record old positions in molecule i
            if(vars->rigid_pairs[i].size()) {
                vars->oldPosition.clear();
                for (auto &a : vars->positionAtom[i]) {
                    std::array<double,3> pos;
                    pos[0]=a[0];
                    pos[1]=a[1];
                    pos[2]=a[2];
                    vars->oldPosition.push_back(pos);
                }
            }

            // update position without constraint
            std::array<double,3> *fi = vars->forceAtom[i].data();
            std::array<double,3> *xi = vars->positionAtom[i].data();
            std::array<double,3> *vi = vars->velocityAtom[i].data();
            int Natom=vars->massAtom[i].size();
            for (int ii=0; ii<Natom; ii++) {
                xi[ii][0] += vi[ii][0] * dt;
                xi[ii][1] += vi[ii][1] * dt;
                xi[ii][2] += vi[ii][2] * dt;
                fi[ii][0]=fi[ii][1]=fi[ii][2]=0.0;
            }

            // calculate constraint force
            if(vars->rigid_pairs[i].size()) shake(i);
        }
    }
    vars->times.tpos+=omp_get_wtime();
}


void
MD::shake(int ID){
    int l=vars->rigid_pairs[ID].size();
    bool converge=false;
    std::array<double,3> *x = vars->positionAtom[ID].data();
    std::array<double,3> *v = vars->velocityAtom[ID].data();
    std::array<double,3> *r = vars->oldPosition.data();
    double *m = vars->massAtom[ID].data();
    while(!converge){
        converge=true;
        for(int k=0;k<l;k++){
            int i=vars->rigid_pairs[ID][k][0];
            int j=vars->rigid_pairs[ID][k][1];
            int ty=vars->rigid_pairs[ID][k][2];
            double dk=vars->btypes[ty].coeff[1];
            double rkx=r[i][0]-r[j][0];
            double rky=r[i][1]-r[j][1];
            double rkz=r[i][2]-r[j][2];
            double rk2=rkx*rkx+rky*rky+rkz*rkz;
            double _rkx=x[i][0]-x[j][0];
            double _rky=x[i][1]-x[j][1];
            double _rkz=x[i][2]-x[j][2];
            double _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dot=rkx*_rkx+rky*_rky+rkz*_rkz;
            double mred_inv=1/m[i]+1/m[j];
            double ramk=(dk*dk-_rk2)/2*mred_inv;
            if(dot>1e-10) ramk/=dot;
            else ramk=0.5*rho_zero+ramk/(rk2*rho_zero)*mred_inv;

            double coeffi=ramk/m[i];
            double coeffj=ramk/m[j];
            x[i][0]+=rkx*coeffi;
            x[i][1]+=rky*coeffi;
            x[i][2]+=rkz*coeffi;
            x[j][0]-=rkx*coeffj;
            x[j][1]-=rky*coeffj;
            x[j][2]-=rkz*coeffj;
            _rkx=x[i][0]-x[j][0];
            _rky=x[i][1]-x[j][1];
            _rkz=x[i][2]-x[j][2];
            _rk2=_rkx*_rkx+_rky*_rky+_rkz*_rkz;
            double dr=_rk2-dk*dk;
            if(dr*dr>tolerance) {
                converge=false;
            }
            if(dr*dr>1e4) converge=false;
        }
    }

}
