#include "../md.hpp"

/////////////////////////////////////////////////////////////////////
/*
	- Randomly arraying gas molecule around an ion with avoiding the
	overlapping. The velocity is picked from  the Maxwell-Boltzumann
	distribution
	- Set ion's center of mass (maybe -> 0), make pair list for initial
	step of simulation, reset the margine size.
*/
/////////////////////////////////////////////////////////////////////

void
MD::initialization_vapor(void) {
	double dis=10;	/*	minimum distance */

    // Set random number generator
	random_device seed;
	mt19937 mt(seed());
	normal_distribution<> distvapor(0.0, sqrt(kb*T/pp->mvapor));
	uniform_real_distribution<double> r(-d_size*0.5,d_size*0.5);

	int NatomOffset=vars->mass.size();

    // Main part, generate random x, y, z positions and calculate minimum gas-gas distance.
	int i=0;
	do {
		std::array<double,3> X;
		std::array<double,3> V;
		std::array<double,3> Force;
		Force[0]=Force[1]=Force[2]=0;
		double min_dis=10000.0;
		X[0]=r(mt), X[1]=r(mt), X[2]=r(mt);
		for(auto &mol : vars->positionAtom){
			for(auto &x : mol){
				double dx=x[0]-X[0];
				double dy=x[1]-X[1];
				double dz=x[2]-X[2];
				adjust_periodic(dx, dy, dz, d_size);
				double d=sqrt(dx*dx+dy*dy+dz*dz);
				if(d<min_dis) min_dis=d; // minimum gas-gas distance
			}
		}
		if(min_dis>dis){
			V[0]=distvapor(mt)*1e-5;
			V[1]=distvapor(mt)*1e-5;
			V[2]=distvapor(mt)*1e-5;
			vars->position.push_back(X);
			vars->velocity.push_back(V);
			vars->force.push_back(Force);
			vars->mass.push_back(pp->Mvapor);
			vars->region.push_back(false);
			vars->group[2].push_back(i+NatomOffset);

			std::vector<std::array<double,3>> dums;
			std::array<double,3> dum;
			dum[0]=dum[1]=dum[2]=0;
			for(auto x : vars->positionVapor) dums.push_back(dum);
			vars->positionAtom.push_back(dums);
			vars->velocityAtom.push_back(dums);
			vars->forceAtom.push_back(dums);
			vars->massAtom.push_back(vars->massVapor);
			vars->typeAtom.push_back(vars->typeVapor);
			vars->chargeAtom.push_back(vars->chargeVapor);
			vars->bonds.push_back(vars->bondsVapor);
			vars->angles.push_back(vars->anglesVapor);
			vars->dihedrals.push_back(vars->dihedralsVapor);
			vars->rigid_pairs.push_back(vars->bondsVaporRigid);

			std::vector<double> dumw;
			for(auto x : vars->positionVapor) dumw.push_back(0);
			vars->weight.push_back(dumw);

			i++;
		}
	} while(i<Nof_around_vapor);
}