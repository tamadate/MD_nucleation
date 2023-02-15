#include "../md.hpp"
/*########################################################################################

-----Periodic conditions-----

#######################################################################################*/


/************************periodec condition for in gas***************************/
void
MD::periodic(void) {
	std::array<double,3> *mol = vars->position.data();
	double HL=d_size*0.5;
	for (auto &a : vars->position) {
		if (a[0] < mol[0][0]-HL) a[0] += d_size;
		if (a[1] < mol[0][1]-HL) a[1] += d_size;
		if (a[2] < mol[0][2]-HL) a[2] += d_size;
		if (a[0] > mol[0][0]+HL) a[0] -= d_size;
		if (a[1] > mol[0][1]+HL) a[1] -= d_size;
		if (a[2] > mol[0][2]+HL) a[2] -= d_size;
	}
}

/************************periodec condition for in gas***************************/
void
MD::boundary_scaling_gas_move(void){

	std::array<double,3> *mol = vars->position.data();
	std::array<double,3> *molv = vars->velocity.data();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;

	for (auto i : vars->gas_out){
		flag=flagx=flagy=flagz=0;
		if (mol[i][0] < preCenter[0]-HL) mol[i][0] += d_size, flagx--, flag++;
		if (mol[i][1] < preCenter[1]-HL) mol[i][1] += d_size, flagy--, flag++;
		if (mol[i][2] < preCenter[2]-HL) mol[i][2] += d_size, flagz--, flag++;
		if (mol[i][0] > preCenter[0]+HL) mol[i][0] -= d_size, flagx++, flag++;
		if (mol[i][1] > preCenter[1]+HL) mol[i][1] -= d_size, flagy++, flag++;
		if (mol[i][2] > preCenter[2]+HL) mol[i][2] -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdist->number>mbdist->vflux.size()*0.9) {mbdist->makeWeightedMB(pp->cgas,pp->mgas,T);}
		    double vx=molv[i][0];
		    double vy=molv[i][1];
		    double vz=molv[i][2];
			double v2 = vx*vx+vy*vy+vz*vz;
			double v = sqrt(v2);
			double vMB = mbdist->vflux[mbdist->number]*1e-5;
			double mod_factor= vMB/v;
		    molv[i][0] *= mod_factor;
		    molv[i][1] *= mod_factor;
		    molv[i][2] *= mod_factor;
		    mbdist->number++;
		}
	}
}

void
MD::boundary_scaling_vapor_move(void){
	std::array<double,3> *mol = vars->position.data();
	std::array<double,3> *molv = vars->velocity.data();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;

	for (auto &i : vars->vapor_out){
		flag=flagx=flagy=flagz=0;
		if (mol[i][0] < preCenter[0]-HL) mol[i][0] += d_size, flagx--, flag++;
		if (mol[i][1] < preCenter[1]-HL) mol[i][1] += d_size, flagy--, flag++;
		if (mol[i][2] < preCenter[2]-HL) mol[i][2] += d_size, flagz--, flag++;
		if (mol[i][0] > preCenter[0]+HL) mol[i][0] -= d_size, flagx++, flag++;
		if (mol[i][1] > preCenter[1]+HL) mol[i][1] -= d_size, flagy++, flag++;
		if (mol[i][2] > preCenter[2]+HL) mol[i][2] -= d_size, flagz++, flag++;
		if (flag>0) {
			if(mbdistV->number>mbdistV->vflux.size()*0.9) {mbdistV->makeWeightedMB(pp->cvapor,pp->mvapor,T);}
		    double vx=molv[i][0];
		    double vy=molv[i][1];
		    double vz=molv[i][2];
			double v2 = vx*vx+vy*vy+vz*vz;
			double v = sqrt(v2);
			double vMB = mbdist->vflux[mbdist->number]*1e-5;
			double mod_factor= vMB/v;
		    molv[i][0] *= mod_factor;
		    molv[i][1] *= mod_factor;
		    molv[i][2] *= mod_factor;
		    mbdistV->number++;
		}
	}
}

void
MD::boundary_scaling_ion_move(void){
	std::array<double,3> *mol = vars->position.data();
	std::array<double,3> *molv = vars->velocity.data();
	double HL=d_size*0.5;
	int flag, flagx, flagy, flagz;
	random_device seed;
	mt19937 mt(seed());
	normal_distribution<> distgas(0.0, sqrt(kb*T/pp->mgas));
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));

	for (auto i : vars->gas_out){
		flag=flagx=flagy=flagz=0;
		if (mol[i][0] < mol[0][0]-HL) mol[i][0] += d_size, flagx--, flag++;
		if (mol[i][1] < mol[0][1]-HL) mol[i][1] += d_size, flagy--, flag++;
		if (mol[i][2] < mol[0][2]-HL) mol[i][2] += d_size, flagz--, flag++;
		if (mol[i][0] > mol[0][0]+HL) mol[i][0] -= d_size, flagx++, flag++;
		if (mol[i][1] > mol[0][1]+HL) mol[i][1] -= d_size, flagy++, flag++;
		if (mol[i][2] > mol[0][2]+HL) mol[i][2] -= d_size, flagz++, flag++;
		if (flag>0) {
			molv[i][0] = distgas(mt) *1e-5;
			molv[i][1] = distgas(mt) *1e-5;
			molv[i][2] = distgas(mt) *1e-5;
		}
	}

	for (auto i : vars->vapor_out){
		flag=flagx=flagy=flagz=0;
		if (mol[i][0] < mol[0][0]-HL) mol[i][0] += d_size, flagx--, flag++;
		if (mol[i][1] < mol[0][1]-HL) mol[i][1] += d_size, flagy--, flag++;
		if (mol[i][2] < mol[0][2]-HL) mol[i][2] += d_size, flagz--, flag++;
		if (mol[i][0] > mol[0][0]+HL) mol[i][0] -= d_size, flagx++, flag++;
		if (mol[i][1] > mol[0][1]+HL) mol[i][1] -= d_size, flagy++, flag++;
		if (mol[i][2] > mol[0][2]+HL) mol[i][2] -= d_size, flagz++, flag++;
		if (flag>0) {
			molv[i][0] = distvapor(mt) *1e-5;
			molv[i][1] = distvapor(mt) *1e-5;
			molv[i][2] = distvapor(mt) *1e-5;
		}
	}

}

void
adjust_periodic(double &dx, double &dy, double &dz, double d_size) {
	const double LH = d_size * 0.5;
	if (dx < -LH)dx += d_size;
	if (dx > LH) dx -= d_size;
	if (dy < -LH)dy += d_size;
	if (dy > LH) dy -= d_size;
	if (dz < -LH)dz += d_size;
	if (dz > LH) dz -= d_size;
}