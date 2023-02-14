#include "potential.hpp"


/////////////////////////////////////////////////////////////////////
/*
	- Calculate force working on gas-gas (LJ)
*/
/////////////////////////////////////////////////////////////////////
void
PotentialGasGas::compute(Variables *vars, FLAG *flags) {
	for(auto &p : pairList){
		int i=p.i;
		int j=p.j;
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *fi = vars->forceAtom[i].data();
		int Ni=vars->positionAtom[i].size();
		std::array<double,3> *xj = vars->positionAtom[j].data();
		std::array<double,3> *fj = vars->forceAtom[j].data();
		int Nj=vars->positionAtom[i].size();
		for (int ii=0; ii<Ni; ii++){
			for (int jj=0; jj<Nj; jj++){
				double dx = xi[ii][0] - xj[jj][0];
				double dy = xi[ii][1] - xj[jj][1];
				double dz = xi[ii][2] - xj[jj][2];
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=vars->typeAtom[i][ii];
				int type2=vars->typeAtom[j][jj];
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double w=vars->weight[i][ii] * vars->weight[0][jj];
				double force_pair = (force_lj)*r2inv*w;
				fi[ii][0] += force_pair * dx;
				fi[ii][1] += force_pair * dy;
				fi[ii][2] += force_pair * dz;
				fj[jj][0] -= force_pair * dx;
				fj[jj][1] -= force_pair * dy;
				fj[jj][2] -= force_pair * dz;
				if(flags->eflag) {
					vars->Utotal.Uvg+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
				}
					//	if(flags->eflag) vars->totalPotential+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					//	vars->totalVirial+=force_lj;
			}
		}
	}
}

void
PotentialGasGas::makePair(Variables *vars) {
    pairList.clear();
	std::array<double,3> *mol=vars->position.data();
	int gs=vars->gas_in.size();
	for (int i=0; i<gs-1; i++){
		for (int j=i+1; j<gs; j++){
			double dx = mol[i][0] - mol[j][0];
			double dy = mol[i][1] - mol[j][1];
			double dz = mol[i][2] - mol[j][2];
			double r2 = (dx * dx + dy * dy + dz * dz);
			if (r2 < ML2){
				Pair p;
				p.i=i;
				p.j=j;
				pairList.push_back(p);
			}
		}
	}

}