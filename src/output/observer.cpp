//------------------------------------------------------------------------
#include "observer.hpp"
//------------------------------------------------------------------------
void
Observer::computeProps(Variables *vars){
	std::array<double,3> *v = vars->velocity.data();

	K_g=0;
	for(auto i : vars->group[1]){
		K_g += v[i][0] * v[i][0] * vars->mass[i];
		K_g += v[i][1] * v[i][1] * vars->mass[i];
		K_g += v[i][2] * v[i][2] * vars->mass[i];
	}
	K_g*= (0.5 * real_to_kcalmol);
	T_g=K_g/double(vars->group[1].size())*coeff;

	K_v=0;
	for(auto i : vars->group[2]){
		K_v += v[i][0] * v[i][0] * vars->mass[i];
		K_v += v[i][1] * v[i][1] * vars->mass[i];
		K_v += v[i][2] * v[i][2] * vars->mass[i];
	}
	K_v*= 0.5 * real_to_kcalmol;
	T_v=K_v/double(vars->group[2].size())*coeff;

	std::array<double,3> *vi = vars->velocityAtom[0].data();
	Kion=0;
	int Nion=vars->velocityAtom[0].size();
	for(int i=0; i<Nion; i++){
		Kion += vi[i][0] * vi[i][0] * vars->massAtom[0][i];
		Kion += vi[i][1] * vi[i][1] * vars->massAtom[0][i];
		Kion += vi[i][2] * vi[i][2] * vars->massAtom[0][i];
	}
	Kion*= 0.5 * real_to_kcalmol;
	Tion=Kion/double(Nion)*coeff;
};

double
Observer::pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T) {
	double phi = 0.0;
	/*const int ps = pairs.size();
	Gas *gases = vars->gases.data();
	for (int k = 0; k < ps; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = gases[j].qx - gases[i].qx;
		double dy = gases[j].qy - gases[i].qy;
		double dz = gases[j].qz - gases[i].qz;
		adjust_periodic(dx, dy, dz);
		double r2 = (dx * dx + dy * dy + dz * dz);
		double r2inv= 1/r2;
		int type1=gases[i].type;
		int type2=gases[j].type;
		if (r2 < CL2){
			double r6inv = r2inv * r2inv * r2inv;
			phi += r6inv * (vars->pair_coeff[type1][type2][0] * r6inv - vars->pair_coeff[type1][type2][1]);
		}
	}*/
	phi = phi * Cpress;
	return  p/T*Treal + phi;
}
