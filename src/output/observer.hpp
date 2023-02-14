#pragma once
#include "../variables.hpp"
//------------------------------------------------------------------------
class Observer {
public:
	const double coeff=kb_real_inv/1.5;
	double K_g;
	double T_g;

	double K_v;
	double T_v;

	double Kion;
	double Tion;

	int Ngas;
	int Nion;
	int Nvap;

	void computeProps(Variables *vars);
	double potential_energy(Variables *vars) {return vars->totalPotential;}
	double pressure(Variables *vars, std::vector<Pair> &pairs, double Treal, double virial,double p,double T);
	double total_energy(Variables *vars, std::vector<Pair> &pairs) {return K_g + potential_energy(vars);}
};
//------------------------------------------------------------------------
