#include "../md.hpp"


void
MD::setDefaultVariables(void){
    dt = 0.5;	/*	fs	*/
	CUTOFF = 20.0;	/*	A	*/
	MARGIN = 10.0;	/*	A	*/
	OBSERVE=10000000;
	T=300;
	p=1e5;
	positionLogStep=0;

	gastype=2;	/*1:He, 2:Ar, 3:N2*/
	step_relax=0;
	Noftimestep=0;
	p=1e5;
	T=300;
	Nof_around_gas=1;
	Nof_around_vapor=1;
	Ecoeff[0]=Ecoeff[1]=Ecoeff[2]=0;
	itime=0;
}