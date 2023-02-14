//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------
void
MD::velocity_scaling(void) {
	obs->computeProps(vars);
	double Tp;
	Tp=300;
//	Tp=2500-2.2e-4*vars->time;
	double ratio=sqrt(Tp/obs->Tion);
	for (auto &a : vars->velocityAtom[0]){
		a[0]*=ratio;
		a[1]*=ratio;
		a[2]*=ratio;
	}
}

void
MD::nosehoover_zeta(void){
	obs->computeProps(vars);
	int g=vars->massAtom[0].size()*3;
	double Q_inv = 0.0001;
	vars->zeta_ion += (obs->Tion - pp->Tnh_ion)*g*kb_real*Q_inv*dt;
}


void
MD::nosehoover_ion(void){
	double Coeff=exp(-vars->zeta_ion*0.5*dt);
	for (auto &a : vars->velocityAtom[0]){
		a[0]*=Coeff;
		a[1]*=Coeff;
		a[2]*=Coeff;
	}
}

void
MD::setNVE(void){
	flags->velocity_scaling=0;
	flags->nose_hoover_ion=0;
	flags->nose_hoover_gas=0;
}


