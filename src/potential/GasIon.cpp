#include "potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on ion-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialGasIon::compute(Variables *vars, FLAG *flags) {
	vars->times.tgi-=omp_get_wtime();
	std::array<double,3> *x = vars->positionAtom[0].data();
	std::array<double,3> *f = vars->forceAtom[0].data();
	for(auto &p : pairList){
		int i=p.i;
		int j=p.j;
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *fi = vars->forceAtom[i].data();
		int Ngas=vars->positionAtom[i].size();
		for (int ii=0; ii<Ngas; ii++){
			double dx = xi[ii][0] - x[j][0];
			double dy = xi[ii][1] - x[j][1];
			double dz = xi[ii][2] - x[j][2];
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			int type1=vars->typeAtom[i][ii];
			int type2=vars->typeAtom[0][j];
			double r6inv = r2inv * r2inv * r2inv;
			double force_pair = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1])*r2inv;
			fi[ii][0] += force_pair * dx;
			fi[ii][1] += force_pair * dy;
			fi[ii][2] += force_pair * dz;
			f[j][0] -= force_pair * dx;
			f[j][1] -= force_pair * dy;
			f[j][2] -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Ugi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
			}
				//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				//	vars->totalVirial+=force_lj;
		}
	}
	vars->times.tgi+=omp_get_wtime();
}

void
PotentialGasIon::makePair(Variables *vars) {
    pairList.clear();
	std::array<double,3> *mol=vars->position.data();
	int Nion=vars->positionAtom[0].size();
	for (auto i : vars->gas_in){
		for(int j=0; j<Nion; j++){
			double dx=mol[i][0]-vars->positionAtom[0][j][0];
			double dy=mol[i][1]-vars->positionAtom[0][j][1];
			double dz=mol[i][2]-vars->positionAtom[0][j][2];
			double r2 = (dx * dx + dy * dy + dz * dz);
			if(r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				pairList.push_back(p);
			}
		}
	}

}