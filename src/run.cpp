//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	- Diffusion coefficient calculation
	- MD simulation performed for (step_relax) steps as a thermal
	relaxation, and main calculation run for (Noftimestep) steps.
	- Among of two calculations, a calculation is re-initialized.
	- Properties (pressure, temperature, energy etc...) are observed
	each (OBSERVE) steps.
*/
/////////////////////////////////////////////////////////////////////
void
MD::run(char** argv) {
/****Thermal relaxation****/
	const int logger=10000;
	for (auto &a : IntraInter) {a->printName();}
	for (auto &a : InterInter) {a->printName();}
	for (itime=0; itime < step_relax; itime++) {
		if (itime%OBSERVE==0) {
			display(1);
			export_dump_close();
			vars->time=(itime*dt);
		}
    	verlet();
		//if ((itime+1)%pp->OBSERVE==0) flags->eflag=1;
	}

	/*
	****Reinitialization****
	- Reset time (time -> 0).
	- Re-initializating the ion and gas positions and velocities.
	Position -> 0 elta, Traslational velocity -> 0
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
	*/
	vars->time=0;
	itime=0;
	for (auto &a : vars->position){
		a[0]-=vars->position[0][0];
		a[1]-=vars->position[0][1];
		a[2]-=vars->position[0][2];
	}
	for (auto &as : vars->positionAtom){
		for (auto &a : as){
			a[0]-=vars->position[0][0];
			a[1]-=vars->position[0][1];
			a[2]-=vars->position[0][2];
		}
	}
	make_pair();
	margin_length = MARGIN;
	setNVE();

	/****Main simulaiton****/
	for (itime=0; itime<Noftimestep; itime++) {
		if(itime%logger==0){
			output();
			vars->time+=(logger*dt);
		}
		if (itime%OBSERVE==0) {
			display(0);
			export_dump_close();
		}
		if ((itime+1)%OBSERVE==0) {
			flags->eflag=1;
			vars->Uzero();
		}
		verlet();
	}
}
