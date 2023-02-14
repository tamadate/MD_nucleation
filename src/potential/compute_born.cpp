#include "potential.hpp"
/*########################################################################################

-----compute intramolecular interaction-----

#######################################################################################*/

/**********************************Force calculation******************************************/
void
PotentialBorn::compute(Variables *vars, FLAG *flags) {
	const int is = vars->massAtom[0].size();
	std::array<double,3> *x = vars->positionAtom[0].data();
	std::array<double,3> *f = vars->forceAtom[0].data();
	vars->times.tion-=omp_get_wtime();
	for(int i=0; i<is-1; i++){
		for(int j=i+1; j<is; j++){
			double dx = x[i][0] - x[j][0];
			double dy = x[i][1] - x[j][1];
			double dz = x[i][2] - x[j][2];
			double rsq = (dx * dx + dy * dy + dz * dz);
			double r2inv = 1/rsq;
			double r6inv = r2inv * r2inv * r2inv;
			double r=sqrt(rsq);
			int type1=vars->typeAtom[0][i];
			int type2=vars->typeAtom[0][j];

			//vars->bornCoeff[type1][type2][0]=A
			//vars->bornCoeff[type1][type2][1]=6C
			//vars->bornCoeff[type1][type2][2]=8D
			//vars->bornCoeff[type1][type2][3]=sigma
			//vars->bornCoeff[type1][type2][4]=1/rho

			double rexp=exp(vars->bornCoeff[type1][type2][4]*(vars->bornCoeff[type1][type2][3]-r));
			double force1 = vars->bornCoeff[type1][type2][4]*vars->bornCoeff[type1][type2][0]*r*rexp;
			double force2 = -vars->bornCoeff[type1][type2][1]*r6inv;
			double force3 = -vars->bornCoeff[type1][type2][2]*r6inv*r2inv;
			double force_coul = qqrd2e*vars->chargeAtom[0][i]*vars->chargeAtom[0][j]*sqrt(r2inv);
			double force_pair = (force1+force2+force3+force_coul)*r2inv;

			f[i][0] += force_pair * dx;
			f[i][1] += force_pair * dy;
			f[i][2] += force_pair * dz;
			f[j][0] -= force_pair * dx;
			f[j][1] -= force_pair * dy;
			f[j][2] -= force_pair * dz;
			if(flags->eflag) {
				vars->Utotal.Uion+=force_coul;
				vars->Utotal.Uion+=rexp*vars->bornCoeff[type1][type2][0];
				vars->Utotal.Uion-=vars->bornCoeff[type1][type2][1]/6.0*r6inv;
				vars->Utotal.Uion-=vars->bornCoeff[type1][type2][2]/8.0*r6inv*r2inv;
			}
			//vars->totalVirial+=force_lj;
		}
	}
	vars->times.tion+=omp_get_wtime();
}
