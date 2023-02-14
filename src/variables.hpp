#pragma once
#include "constants.hpp"
#include "flags.hpp"

//------------------------------------------------------------------------
class Variables {
public:

	Variables(void);
	~Variables(void){
		time=0;
		Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
		flags = new FLAG();
	};

	FLAG *flags;

	/*variables*/
  	int Nth;
	std::vector<std::array<double,3>> position;
	std::vector<std::array<double,3>> velocity;
	std::vector<std::array<double,3>> force;
	std::vector<double> mass;
	std::vector<double> type;
	std::vector<bool> region;

	std::vector<std::vector<std::array<double,3>>> positionAtom;
	std::vector<std::vector<std::array<double,3>>> velocityAtom;
	std::vector<std::vector<std::array<double,3>>> forceAtom;
	std::vector<std::vector<double>> massAtom;
	std::vector<std::vector<double>> chargeAtom;
	std::vector<std::vector<int>> typeAtom;
	std::vector<std::vector<double>> weight;
	std::vector<std::vector<std::array<double,3>>> bonds;
	std::vector<std::vector<std::array<double,4>>> angles;
	std::vector<std::vector<std::array<double,5>>> dihedrals;
	std::vector<std::vector<std::array<double,3>>> rigid_pairs;		// SHAKE
    std::vector<std::array<double,3>> oldPosition;					// For a molecule

	std::vector<std::array<double,3>> positionVapor;
	std::vector<double> massVapor;
	std::vector<double> chargeVapor;
	std::vector<int> typeVapor;
	std::vector<std::array<double,3>> bondsVapor;
	std::vector<std::array<double,4>> anglesVapor;
	std::vector<std::array<double,5>> dihedralsVapor;
	std::vector<std::array<double,3>> bondsVaporRigid;		// SHAKE

	std::vector<std::array<double,3>> positionGas;
	std::vector<std::array<double,3>> velocityGas;
	std::vector<std::array<double,3>> forceGas;
	std::vector<int> typeGas;
	std::vector<double> massGas;
	std::vector<double> chargeGas;
	std::vector<std::array<double,3>> bondsGas;
	std::vector<std::array<double,3>> bondsGasRigid;		// SHAKE

	std::array<std::vector<int>,3> group;
	void setCrossPotentials(void);

	double time;
	double zeta_ion;
	double zeta_gas;

	Potentials Utotal;
	Times times;

	void Uzero(void)	{
		Utotal.Uion=Utotal.Ugas=Utotal.Uvap=Utotal.Ugi=Utotal.Ugg=Utotal.Uvg=Utotal.Uvi=Utotal.Uvv=0;
  	}

  	void tzero(void)	{times.tion=times.tgas=times.tvap=times.tgi=times.tvv=times.tvg=times.tvi=times.tpair=0;}
	double Usum(void)	{return Utotal.Uion+Utotal.Ugas+Utotal.Uvap+Utotal.Ugi+Utotal.Ugg+Utotal.Uvg+Utotal.Uvi+Utotal.Uvv;}

	std::vector<int> gas_in;	/*	gas list around ion1	*/
	std::vector<int> gas_out;	/*	gas list far from ion1	*/
	std::vector<int> vapor_in;	/*	gas list around ion1	*/
	std::vector<int> vapor_out;	/*	gas list far from ion1	*/

	std::vector<Atom_type> atypes;
	std::vector<Bond_type> btypes;
	std::vector<Angle_type> ctypes;
	std::vector<Dihedral_type> dtypes;

	std::vector<vector<vector<double>>> pair_coeff;
	double bornCoeff[2][2][5];

	void setGasPotentials(void);
	void setBMHPotential(void);
	void setCrossPotentials(int Nion,int Nvapor);

	/*initialization and export to dump file*/
	void read_initial(char* ionFile, char* vaporFile);
	void readIonFile(char* infile);
	void readVaporFile(char* infile);
	void ionInitialVelocity(double T);
	void ionRotation(void);
	double totalPotential;
	double totalVirial;

	void ROTATION(double &X, double &Y, double &Z, double A, double B, double C, double x, double y, double z){
		X = cos(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+sin(C)*(y*cos(A)+z*sin(A));
		Y = -sin(C)*(x*cos(B)+y*sin(A)*sin(B)-z*cos(A)*sin(B))+cos(C)*(y*cos(A)+z*sin(A));
		Z = x*sin(B)-y*sin(A)*cos(B)+z*cos(A)*cos(B);
	}

private:
};
//------------------------------------------------------------------------
