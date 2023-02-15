#include "potential.hpp"


void
PotentialVaporIon::compute(Variables *vars, FLAG *flags) {
	vars->times.tvi-=omp_get_wtime();
	std::array<double,3> *xj = vars->positionAtom[0].data();
	std::array<double,3> *fj = vars->forceAtom[0].data();
	int Nj=vars->positionAtom[0].size();
	for(auto i : vars->vapor_in){
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *fi = vars->forceAtom[i].data();
		int Ni = vars->massAtom[i].size();
		for(int ii=0;ii<Ni;ii++){
			for(int jj=0;jj<Nj;jj++){
				double dx = xi[ii][0] - xj[jj][0];
				double dy = xi[ii][1] - xj[jj][1];
				double dz = xi[ii][2] - xj[jj][2];
				double rsq = (dx * dx + dy * dy + dz * dz);
				double r2inv = 1/rsq;
				int type1=vars->typeAtom[i][ii];
				int type2=vars->typeAtom[0][jj];
				double r6inv = r2inv * r2inv * r2inv;
				double force_lj = r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
				double force_coul = qqrd2e * vars->chargeAtom[i][ii] * vars->chargeAtom[0][jj] * sqrt(r2inv);
				double w=vars->weight[i][ii] * vars->weight[0][jj];
				double force_pair = (force_lj + force_coul)*r2inv*w;
				fi[ii][0] += force_pair * dx;
				fi[ii][1] += force_pair * dy;
				fi[ii][2] += force_pair * dz;
				fj[jj][0] -= force_pair * dx;
				fj[jj][1] -= force_pair * dy;
				fj[jj][2] -= force_pair * dz;
				if(flags->eflag) {
					vars->Utotal.Uvi+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					vars->Utotal.Uvi+=force_coul;
				}
				//	vars->totalVirial+=force_lj;
			}
		}
	}
	vars->times.tvi+=omp_get_wtime();
}
