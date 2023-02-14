//------------------------------------------------------------------------
#include "../md.hpp"
//------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////////
/*
	make pair list
	- always gas-ion pair list was calculated
	- if gas-gas interaction is ON, make pair_gasgas
	- get max velocity of gas molecule vmax
	- set update loop = margine_length/(vmax*dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::make_pair(void){
	vars->times.tpair-=omp_get_wtime();
	std::array<double,3> *mol=vars->position.data();
	std::array<double,3> *molv=vars->velocity.data();
	for (auto i : vars->gas_out) {
		mol[i][0] += molv[i][0]*dt*loopPair;
		mol[i][1] += molv[i][1]*dt*loopPair;
		mol[i][2] += molv[i][2]*dt*loopPair;
	}
	for (auto i : vars->vapor_out) {
		mol[i][0] += molv[i][0]*dt*loopPair;
		mol[i][1] += molv[i][1]*dt*loopPair;
		mol[i][2] += molv[i][2]*dt*loopPair;
	}
	/*boundary_scaling_gas_move();
	boundary_scaling_vapor_move();
	boundary_scaling_ion_move();
	pre_ion[0]=mol[0][0];
	pre_ion[1]=mol[0][1];
	pre_ion[2]=mol[0][2];*/

	updateInCenters();
	update_gas_in();
	update_vapor_in();
	
	for (auto &a : InterInter) a->makePair(vars);

	//	set number of steps to update pair list
	double vmax2 = 0.0;
	for (auto &v : vars->velocity) {
		double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
		if (vmax2 < v2) vmax2 = v2;
	}
	double vmax = sqrt(vmax2);

	loop_update=margin_length/(vmax*dt);
	loopPair=0;

	vars->times.tpair+=omp_get_wtime();
}


/////////////////////////////////////////////////////////////////////
/*
	check necessity of the pair list updating
*/
/////////////////////////////////////////////////////////////////////
void
MD::check_pairlist(void){
	loopPair++;
	if(loopPair>loop_update) make_pair();
//	if(flags->force_sw==1) sw->check_pairlist(vars);
//	if(flags->force_ters==1) ters->check_pairlist(vars);
}
