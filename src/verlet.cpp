//------------------------------------------------------------------------
#include "md.hpp"
//------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////////
/*
	- Compute a domain (NVE or NVT)
	- v(t) -> v(t+dt/2)
	- Calculate velocity and potition of ion's center of mass.
	- Temperature control (velocity scaling, this case is NVT)
	- r(t) -> r(t+dt)
	- apply periodic boundary condition
	- Determine whethere update pair list or not
	- Calculate ion intraatmic interaction
	- Calculate ion-gas interatominc interaction
	- v(t+dt/2) -> v(t)
*/
/////////////////////////////////////////////////////////////////////
void
MD::verlet(void) {
	velocity_calculation(); //	v(t) -> v(t+dt/2) using F(x(t))
	if(flags->velocity_scaling==1)	velocity_scaling();

	// updatePositionSHAKE();
	update_position();

	if(flags->nose_hoover_ion==1)	nosehoover_zeta();
	//if(flags->nose_hoover_gas==1)	nosehoover_zeta_gas();
	check_pairlist();
	vars->totalVirial=0;
	for (auto &a : InterInter) a->compute(vars,flags);
	for (auto &a : IntraInter) a->compute(vars,flags);

	velocity_calculation();	//	v(t+dt/2) -> v(t+dt) using F(x(t+dt))

	if(flags->nose_hoover_ion==1)	nosehoover_ion();
	//(flags->nose_hoover_gas==1)	nosehoover_gas();
	flags->eflag=0;
}

/////////////////////////////////////////////////////////////////////
/*
	- Update velocity (half of a time step, dt/2)
*/
/////////////////////////////////////////////////////////////////////
void
MD::velocity_calculation(void) {
	vars->times.tvel-=omp_get_wtime();
	double const Coeff=0.5*dt*4.184e-4;
	int Nmol=vars->mass.size();
	for (int i=0; i<Nmol; i++){
		if(vars->region[i]){
			int Natom=vars->massAtom[i].size();
			std::array<double,3> *vi = vars->velocityAtom[i].data();
			std::array<double,3> *fi = vars->forceAtom[i].data();
			for (int ii=0; ii<Natom; ii++){
				double Coeff2=Coeff/vars->massAtom[i][ii];
				vi[ii][0] += fi[ii][0] *Coeff2;
				vi[ii][1] += fi[ii][1] *Coeff2;
				vi[ii][2] += fi[ii][2] *Coeff2;
			}
		}
	}
	vars->times.tvel+=omp_get_wtime();
}

/////////////////////////////////////////////////////////////////////
/*
	- Update velocity (a time step, dt)
*/
/////////////////////////////////////////////////////////////////////
void
MD::update_position(void) {
	vars->times.tpos-=omp_get_wtime();
	int Nmol=vars->mass.size();
	for (int i=0; i<Nmol; i++){
		if(vars->region[i]){
			int Natom=vars->massAtom[i].size();
			std::array<double,3> *xi = vars->positionAtom[i].data();
			std::array<double,3> *vi = vars->velocityAtom[i].data();
			std::array<double,3> *fi = vars->forceAtom[i].data();
			for (int ii=0; ii<Natom; ii++){
				xi[ii][0] += vi[ii][0] * dt;
				xi[ii][1] += vi[ii][1] * dt;
				xi[ii][2] += vi[ii][2] * dt;
				fi[ii][0]=fi[ii][1]=fi[ii][2]=0.0;
			}
		}
	}
	int is=vars->positionAtom[0].size();
	std::array<double,3> &x = vars->position[0];
	std::array<double,3> &v = vars->velocity[0];
	std::array<double,3> *xi = vars->positionAtom[0].data();
	std::array<double,3> *vi = vars->velocityAtom[0].data();
	x[0]=0;	x[1]=0;	x[2]=0;
	v[0]=0;	v[1]=0;	v[2]=0;
	for (int ii=0; ii<is; ii++){
		x[0]+=xi[ii][0]*vars->massAtom[0][ii];
		x[1]+=xi[ii][1]*vars->massAtom[0][ii];
		x[2]+=xi[ii][2]*vars->massAtom[0][ii];
		v[0]+=vi[ii][0]*vars->massAtom[0][ii];
		v[1]+=vi[ii][1]*vars->massAtom[0][ii];
		v[2]+=vi[ii][2]*vars->massAtom[0][ii];
	}
	double mInv=1/vars->mass[0];
	x[0]*=mInv;	x[1]*=mInv;	x[2]*=mInv;
	v[0]*=mInv;	v[1]*=mInv;	v[2]*=mInv;

	for (int I=1; I<3; I++){
		for (auto i : vars->group[I]){
			if(vars->region[i]){
				int Natom=vars->massAtom[i].size();
				std::array<double,3> *xi = vars->positionAtom[i].data();
				std::array<double,3> *vi = vars->velocityAtom[i].data();
				std::array<double,3> *fi = vars->forceAtom[i].data();
				std::array<double,3> &x0 = vars->position[0];
				for (int ii=0; ii<Natom; ii++){
					double dx = xi[ii][0] - x0[0];
					double dy = xi[ii][1] - x0[1];
					double dz = xi[ii][2] - x0[2];
					double r2=dx*dx+dy*dy+dz*dz;
					if(r2<Ri2) {
						vars->weight[i][ii]=1;
					} else if(r2>RI2) {
						vars->weight[i][ii]=0;
					} else {
						double dr=(sqrt(r2)-Rcenter);
						double C=cos(M_PI/Overlap*0.5*dr);
						vars->weight[i][ii]=C*C;
					}
				}
			}
		}
	}
	vars->times.tpos+=omp_get_wtime();
}


