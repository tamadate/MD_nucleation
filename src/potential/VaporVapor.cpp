#include "potential.hpp"


void
PotentialVaporVapor::compute(Variables *vars, FLAG *flags) {
	vars->times.tvv-=omp_get_wtime();
	for(auto &p : pairList){
		int i=p.i;
		int j=p.j;
		std::array<double,3> *xi = vars->positionAtom[i].data();
		std::array<double,3> *fi = vars->forceAtom[i].data();
		int Ni=vars->positionAtom[i].size();
		std::array<double,3> *xj = vars->positionAtom[j].data();
		std::array<double,3> *fj = vars->forceAtom[j].data();
		int Nj=vars->positionAtom[j].size();
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
				double force_coul = qqrd2e * vars->chargeAtom[i][ii] * vars->chargeAtom[j][jj] * sqrt(r2inv);
				double w=vars->weight[i][ii] * vars->weight[j][jj];
				double force_pair = (force_lj + force_coul)*r2inv*w;
				fi[ii][0] += force_pair * dx;
				fi[ii][1] += force_pair * dy;
				fi[ii][2] += force_pair * dz;
				fj[jj][0] -= force_pair * dx;
				fj[jj][1] -= force_pair * dy;
				fj[jj][2] -= force_pair * dz;
				if(flags->eflag) {
					vars->Utotal.Uvv+=r6inv * (vars->pair_coeff[type1][type2][0]/12.0 * r6inv - vars->pair_coeff[type1][type2][1]/6.0);
					vars->Utotal.Uvv+=force_coul;
				}
				//	vars->totalVirial+=force_lj;
			}
		}
	}
	vars->times.tvv+=omp_get_wtime();
}

void
PotentialVaporVapor::makePair(Variables *vars) {
    pairList.clear();
	const int vs = vars->vapor_in.size();
	if(vs>1){
		Pair p;
		for(int i1=0; i1<vs-1; i1++){
			p.i=vars->vapor_in[i1];
			for(int i2=i1+1; i2<vs; i2++){
				p.j=vars->vapor_in[i2];
				pairList.push_back(p);
			}
		}
	}
}